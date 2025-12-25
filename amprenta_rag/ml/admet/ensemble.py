from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from amprenta_rag.logging_utils import get_logger


logger = get_logger(__name__)


def _require_xgboost():
    try:
        from xgboost import XGBClassifier, XGBRegressor  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("xgboost is required for BootstrapEnsemble (pip install xgboost)") from e
    return XGBClassifier, XGBRegressor


@dataclass(frozen=True)
class EnsembleArtifact:
    models: List[Any]
    training_centroid: np.ndarray
    n_models: int
    feature_dim: int
    task: str  # "classification" | "regression"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "models": self.models,
            "training_centroid": self.training_centroid,
            "n_models": int(self.n_models),
            "feature_dim": int(self.feature_dim),
            "task": self.task,
        }


class BootstrapEnsemble:
    """Bootstrap ensemble over XGBoost models.

    - For classification (binary y in {0,1}), `predict_proba()` returns mean/std of class-1 probability.
    - For regression, `predict()` returns mean/std of predicted value.

    The `training_centroid` is stored for applicability checks (mean fingerprint vector).
    """

    def __init__(self, n_models: int = 5, base_params: Optional[dict] = None):
        self.n_models = int(n_models)
        self.base_params = dict(base_params or {})
        self.models: List[Any] = []
        self.training_centroid: Optional[np.ndarray] = None
        self.feature_dim: Optional[int] = None
        self.task: Optional[str] = None

        # test/debug hooks
        self.bootstrap_indices_: List[np.ndarray] = []

    def fit(self, X: np.ndarray, y: np.ndarray) -> "BootstrapEnsemble":
        if self.n_models <= 0:
            raise ValueError("n_models must be > 0")
        if X is None or y is None:
            raise ValueError("X and y are required")
        if len(X) != len(y):
            raise ValueError("X and y must have the same number of rows")
        if len(X) < 10:
            raise ValueError("Need at least 10 samples for bootstrap training")

        X = np.asarray(X)
        y = np.asarray(y)
        self.feature_dim = int(X.shape[1])
        self.training_centroid = X.mean(axis=0).astype(np.float32, copy=False)

        # Determine task
        uniq = np.unique(y[~np.isnan(y)] if np.issubdtype(y.dtype, np.floating) else y)
        is_binary = set(uniq.tolist()).issubset({0, 1}) and len(uniq) <= 2
        self.task = "classification" if is_binary else "regression"

        XGBClassifier, XGBRegressor = _require_xgboost()

        rng = np.random.default_rng(42)
        self.models = []
        self.bootstrap_indices_ = []

        for i in range(self.n_models):
            idx = rng.integers(0, len(X), size=len(X), endpoint=False)
            self.bootstrap_indices_.append(idx)
            Xb = X[idx]
            yb = y[idx]

            params = dict(self.base_params)
            # keep training fast by default unless caller overrides
            params.setdefault("n_estimators", 200)
            params.setdefault("max_depth", 6)
            params.setdefault("learning_rate", 0.05)
            params.setdefault("subsample", 0.9)
            params.setdefault("colsample_bytree", 0.9)
            params.setdefault("random_state", 1000 + i)
            params.setdefault("n_jobs", 1)

            if self.task == "classification":
                params.setdefault("eval_metric", "logloss")
                model = XGBClassifier(**params)
                model.fit(Xb, yb)
            else:
                model = XGBRegressor(**params)
                model.fit(Xb, yb)

            self.models.append(model)

        # correlation check on training predictions
        try:
            preds = self._pred_matrix(X)
            self._warn_if_high_corr(preds, threshold=0.95)
        except Exception as e:  # noqa: BLE001
            logger.debug("Correlation check skipped/failed: %r", e)

        return self

    def _pred_matrix(self, X: np.ndarray) -> np.ndarray:
        if not self.models:
            raise ValueError("Ensemble not fit")
        if self.task == "classification":
            probs = []
            for m in self.models:
                p = m.predict_proba(X)[:, 1]
                probs.append(np.asarray(p, dtype=float))
            return np.vstack(probs)
        preds = []
        for m in self.models:
            p = m.predict(X)
            preds.append(np.asarray(p, dtype=float))
        return np.vstack(preds)

    @staticmethod
    def _warn_if_high_corr(preds: np.ndarray, threshold: float = 0.95) -> None:
        """Warn if any pair of models is overly correlated."""
        preds = np.asarray(preds, dtype=float)
        if preds.ndim != 2 or preds.shape[0] < 2 or preds.shape[1] < 2:
            return
        corr = np.corrcoef(preds)
        # ignore diagonal
        n = corr.shape[0]
        for i in range(n):
            corr[i, i] = 0.0
        max_corr = float(np.nanmax(np.abs(corr)))
        if max_corr > float(threshold):
            logger.warning("High pairwise ensemble correlation detected: max_corr=%.3f", max_corr)

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if not self.models or not self.task:
            raise ValueError("Ensemble not fit")
        X = np.asarray(X)

        if self.task == "classification":
            return self.predict_proba(X)

        preds = np.vstack([np.asarray(m.predict(X), dtype=float) for m in self.models])
        return preds.mean(axis=0), preds.std(axis=0)

    def predict_proba(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if not self.models or not self.task:
            raise ValueError("Ensemble not fit")
        if self.task != "classification":
            raise ValueError("predict_proba() only valid for classification ensembles")
        X = np.asarray(X)
        probs = np.vstack([np.asarray(m.predict_proba(X)[:, 1], dtype=float) for m in self.models])
        mean = probs.mean(axis=0)
        std = probs.std(axis=0)
        # guard numerical noise
        mean = np.clip(mean, 0.0, 1.0)
        return mean, std

    def to_artifact(self) -> Dict[str, Any]:
        if self.feature_dim is None or self.training_centroid is None or self.task is None:
            raise ValueError("Ensemble not fit")
        art = EnsembleArtifact(
            models=self.models,
            training_centroid=self.training_centroid,
            n_models=self.n_models,
            feature_dim=self.feature_dim,
            task=self.task,
        )
        return art.to_dict()


