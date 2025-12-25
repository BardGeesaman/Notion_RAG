from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


def _ensure_repo_on_path() -> None:
    # When executed as a script ("python scripts/xyz.py"), Python's sys.path[0]
    # is the scripts/ directory, not the repo root. Ensure repo root is importable.
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


def _require_numpy():
    import numpy as np  # type: ignore

    return np


def _require_joblib():
    try:
        import joblib  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("joblib is required (pip install joblib)") from e
    return joblib


def _require_sklearn_metrics_and_cv():
    try:
        from sklearn.metrics import (  # type: ignore
            brier_score_loss,
            mean_squared_error,
            r2_score,
            roc_auc_score,
        )
        from sklearn.model_selection import RandomizedSearchCV  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("scikit-learn is required (pip install scikit-learn)") from e
    return roc_auc_score, brier_score_loss, mean_squared_error, r2_score, RandomizedSearchCV


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _coverage_90(y_true, mean, std) -> float:
    import numpy as np

    z = 1.645  # approx 90% normal interval
    lo = mean - z * std
    hi = mean + z * std
    y_true = np.asarray(y_true, dtype=float)
    return float(np.mean((y_true >= lo) & (y_true <= hi)))


def _print_cfg(cfg: Dict[str, Any]) -> None:
    print(json.dumps(cfg, indent=2, sort_keys=True))


def main() -> None:
    _ensure_repo_on_path()

    ap = argparse.ArgumentParser(description="Train ADMET bootstrap ensembles from curated ChEMBL parquet splits.")
    ap.add_argument("--data-dir", default="data/admet", help="Directory with curated parquet splits")
    ap.add_argument("--endpoints", default="herg,logs,logp", help="Comma-separated endpoints")
    ap.add_argument("--n-models", type=int, default=5, help="Ensemble size")
    ap.add_argument("--output-dir", default="models/admet", help="Output directory for artifacts")
    ap.add_argument("--tune", action="store_true", help="Enable hyperparameter tuning (slow)")
    ap.add_argument("--dry-run", action="store_true", help="Print config and exit without training")
    args = ap.parse_args()

    endpoints = [e.strip().lower() for e in str(args.endpoints).split(",") if e.strip()]
    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)

    cfg = {
        "trained_at": _now_iso(),
        "data_dir": str(data_dir),
        "output_dir": str(out_dir),
        "endpoints": endpoints,
        "n_models": int(args.n_models),
        "tune": bool(args.tune),
        "dry_run": bool(args.dry_run),
    }

    if args.dry_run:
        # No heavy imports here; just show what would run.
        from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader

        loader = ChEMBLDatasetLoader(str(data_dir))
        available = loader.list_endpoints()
        cfg["available_endpoints"] = available
        _print_cfg(cfg)
        return

    # Heavy imports and actual training below.
    np = _require_numpy()
    joblib = _require_joblib()
    roc_auc_score, brier_score_loss, mean_squared_error, r2_score, RandomizedSearchCV = _require_sklearn_metrics_and_cv()

    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper, reliability_diagram
    from amprenta_rag.ml.admet.datasets import ChEMBLDatasetLoader
    from amprenta_rag.ml.admet.ensemble import BootstrapEnsemble, _require_xgboost
    from amprenta_rag.ml.registry import get_registry

    out_dir.mkdir(parents=True, exist_ok=True)
    loader = ChEMBLDatasetLoader(str(data_dir))
    registry = get_registry()

    XGBClassifier, XGBRegressor = _require_xgboost()

    for endpoint in endpoints:
        X_train, y_train, _ = loader.load_endpoint(endpoint, "train")
        X_cal, y_cal, _ = loader.load_endpoint(endpoint, "cal")
        X_test, y_test, _ = loader.load_endpoint(endpoint, "test")

        if len(X_train) == 0 or len(X_test) == 0:
            raise ValueError(f"Empty dataset for endpoint={endpoint}; ensure parquet files exist and SMILES are valid")

        is_classification = set(np.unique(y_train).tolist()).issubset({0.0, 1.0})
        model_type = "admet_classification" if is_classification else "admet_regression"

        base_params: Dict[str, Any] = {}
        if args.tune:
            # Tune a single base model, then apply best params to ensemble.
            if is_classification:
                base = XGBClassifier(eval_metric="logloss", n_jobs=1)
            else:
                base = XGBRegressor(n_jobs=1)

            space = {
                "max_depth": [3, 5, 7, 10],
                "learning_rate": [0.01, 0.05, 0.1, 0.2],
                "n_estimators": [50, 100, 200, 500],
            }
            search = RandomizedSearchCV(
                base,
                param_distributions=space,
                n_iter=20,
                cv=3,
                random_state=42,
                n_jobs=1,
            )
            search.fit(X_train, y_train)
            base_params = dict(search.best_params_ or {})

        ens = BootstrapEnsemble(n_models=int(args.n_models), base_params=base_params).fit(X_train, y_train)

        # Predictions + uncertainty
        if is_classification:
            p_cal, p_cal_std = ens.predict_proba(X_cal)
            p_test, p_test_std = ens.predict_proba(X_test)
            cal = CalibrationWrapper(method="isotonic").fit(p_cal, y_cal)
            p_test_cal = cal.calibrate(p_test)

            auc = float(roc_auc_score(y_test, p_test_cal))
            brier = float(brier_score_loss(y_test, p_test_cal))
            ece = float(reliability_diagram(p_test_cal, y_test, n_bins=10)["ece"])
            coverage = float(_coverage_90(y_test, p_test_cal, p_test_std))

            metrics = {"auc": auc, "brier": brier, "ece": ece, "coverage_90": coverage}
            calibrator_obj: Optional[Any] = cal
        else:
            yhat_cal, yhat_cal_std = ens.predict(X_cal)
            _ = yhat_cal_std  # currently unused; placeholder for regression calibration ideas
            calibrator_obj = None

            yhat, yhat_std = ens.predict(X_test)
            rmse = float(np.sqrt(mean_squared_error(y_test, yhat)))
            r2 = float(r2_score(y_test, yhat))
            coverage = float(_coverage_90(y_test, yhat, yhat_std))

            metrics = {"rmse": rmse, "r2": r2, "coverage_90": coverage}

        app = ApplicabilityChecker(threshold=0.3).fit(X_train)

        artifact = {
            "ensemble": ens.to_artifact(),
            "calibrator": calibrator_obj,
            "applicability": {"centroid": app.training_centroid, "threshold": float(app.threshold)},
            "feature_schema": {"morgan_bits": 2048, "rdkit_descs": 6},
            "metrics": metrics,
            "endpoint": endpoint,
            "trained_at": _now_iso(),
        }

        artifact_path = out_dir / f"{endpoint}_ensemble.joblib"
        joblib.dump(artifact, artifact_path)

        registry.register_model(
            name=f"admet_{endpoint}_ensemble",
            version="1.0.0",
            model_type=model_type,
            framework="xgboost_ensemble",
            model_object=artifact,
            metrics=metrics,
            description=f"Bootstrap ensemble + calibration + applicability for {endpoint}",
        )

        print(f"[OK] Trained + saved {endpoint}: {artifact_path}")


if __name__ == "__main__":
    main()


