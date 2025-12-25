"""Per-target QSAR training using the ADMET ensemble infrastructure."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ml.admet.applicability import ApplicabilityChecker
from amprenta_rag.ml.admet.calibration import CalibrationWrapper
from amprenta_rag.ml.admet.ensemble import BootstrapEnsemble
from amprenta_rag.ml.qsar.datasets import TargetDatasetLoader
from amprenta_rag.ml.registry import get_registry


logger = get_logger(__name__)


def _require_sklearn_split_and_metrics():
    try:
        from sklearn.metrics import accuracy_score, roc_auc_score  # type: ignore
        from sklearn.model_selection import train_test_split  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("scikit-learn is required for QSAR training (pip install scikit-learn)") from e
    return train_test_split, roc_auc_score, accuracy_score


@dataclass(frozen=True)
class QSARTrainingResult:
    artifact: Dict[str, Any]
    metrics: Dict[str, float]
    model_name: str
    version: str
    metadata: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model_name": self.model_name,
            "version": self.version,
            "metrics": self.metrics,
            "metadata": self.metadata,
            "artifact": self.artifact,
        }


def train_target_model(
    target: str,
    source: str = "chembl",
    threshold_nm: float = 1000,
    n_models: int = 5,
    calibration_method: str = "isotonic",
    register: bool = True,
) -> Dict[str, Any]:
    """Train a per-target QSAR binary classifier (IC50 thresholding).

    Returns metrics + model info. Optionally registers artifact in MLModelRegistry.
    """
    train_test_split, roc_auc_score, accuracy_score = _require_sklearn_split_and_metrics()

    loader = TargetDatasetLoader()
    logger.info("[QSAR] Loading dataset target=%s source=%s threshold_nm=%.1f", target, source, float(threshold_nm))
    X, y, smiles, ds_meta = loader.load_target(
        target=target,
        source=source,
        threshold_nm=float(threshold_nm),
        min_compounds=100,
        min_active_ratio=0.2,
    )
    X = np.asarray(X, dtype=np.float32)
    y = np.asarray(y, dtype=np.int64).reshape(-1)
    if X.ndim != 2 or y.ndim != 1 or X.shape[0] != y.shape[0]:
        raise ValueError("Invalid dataset shapes for QSAR training")

    n = int(X.shape[0])
    pos = int((y == 1).sum())
    neg = int((y == 0).sum())
    logger.info("[QSAR] Loaded %d samples (pos=%d neg=%d active_ratio=%.3f)", n, pos, neg, float(pos / max(n, 1)))

    # Split 70/15/15 with stratification if possible.
    strat = y if (pos > 0 and neg > 0) else None
    X_train, X_tmp, y_train, y_tmp = train_test_split(
        X, y, test_size=0.30, random_state=42, stratify=strat
    )
    strat_tmp = y_tmp if (int((y_tmp == 1).sum()) > 0 and int((y_tmp == 0).sum()) > 0) else None
    X_cal, X_test, y_cal, y_test = train_test_split(
        X_tmp, y_tmp, test_size=0.50, random_state=42, stratify=strat_tmp
    )

    # Class imbalance scaling (avoid div-by-zero).
    n_pos = int((y_train == 1).sum())
    n_neg = int((y_train == 0).sum())
    scale_pos_weight = float(n_neg / max(n_pos, 1))
    base_params = {"scale_pos_weight": scale_pos_weight}

    logger.info(
        "[QSAR] Training ensemble n_models=%d scale_pos_weight=%.3f train=%d cal=%d test=%d",
        int(n_models),
        float(scale_pos_weight),
        int(X_train.shape[0]),
        int(X_cal.shape[0]),
        int(X_test.shape[0]),
    )

    ensemble = BootstrapEnsemble(n_models=int(n_models), base_params=base_params).fit(X_train, y_train)

    # Calibrate on calibration split (uses ensemble mean prob as input).
    cal_mean, _cal_std = ensemble.predict_proba(X_cal)
    calibrator = CalibrationWrapper(method=str(calibration_method)).fit(cal_mean, y_cal)

    # Applicability checker on train.
    app = ApplicabilityChecker(threshold=0.3).fit(X_train)

    # Evaluate on test.
    test_mean, test_std = ensemble.predict_proba(X_test)
    test_cal = calibrator.calibrate(test_mean)

    auc = float(roc_auc_score(y_test, test_cal)) if len(np.unique(y_test)) > 1 else float("nan")
    acc = float(accuracy_score(y_test, (test_cal >= 0.5).astype(int)))
    ece = float(calibrator.compute_ece(test_cal, y_test, n_bins=10))

    metrics: Dict[str, float] = {
        "auc": auc,
        "accuracy": acc,
        "ece": ece,
    }
    logger.info("[QSAR] Done target=%s auc=%.3f acc=%.3f ece=%.3f", target, auc, acc, ece)

    artifact: Dict[str, Any] = {
        "ensemble": ensemble.to_artifact(),
        "calibrator": calibrator,
        "applicability": {"threshold": float(app.threshold), "centroid": app.training_centroid},
        "metadata": {
            "target": target,
            "threshold_nm": float(threshold_nm),
            "train_size": int(X_train.shape[0]),
            "calibration_size": int(X_cal.shape[0]),
            "test_size": int(X_test.shape[0]),
            "test_auc": float(auc),
            "source": source,
            "scale_pos_weight": float(scale_pos_weight),
            "dataset_meta": ds_meta,
        },
    }

    model_name = f"qsar_{target}_ensemble"
    version = "1.0.0"
    if register:
        reg = get_registry()
        try:
            reg.register_model(
                name=model_name,
                version=version,
                model_type="qsar_classification",
                framework="xgboost_ensemble",
                model_object=artifact,
                metrics=metrics,
                description=f"Per-target QSAR ensemble for {target} (IC50<{float(threshold_nm)} nM).",
            )
        except Exception as e:  # noqa: BLE001
            logger.warning("[QSAR] Registry registration failed for %s: %s", model_name, e)
            # Keep training result usable even if registry unavailable.
            artifact["metadata"]["registration_error"] = str(e)

    return QSARTrainingResult(
        artifact=artifact,
        metrics=metrics,
        model_name=model_name,
        version=version,
        metadata=artifact["metadata"],
    ).to_dict()


__all__ = ["train_target_model"]


