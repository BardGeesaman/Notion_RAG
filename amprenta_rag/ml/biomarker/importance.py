"""Cross-validated feature importance ranking for biomarker discovery."""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import KFold

    SKLEARN_AVAILABLE = True
except Exception:  # noqa: BLE001
    SKLEARN_AVAILABLE = False


class CVFeatureImportance:
    """KFold CV + RandomForestClassifier feature importances."""

    def __init__(self, n_folds: int = 5, n_estimators: int = 100):
        if not SKLEARN_AVAILABLE:
            raise ImportError(
                "scikit-learn is required for CVFeatureImportance (RandomForestClassifier). "
                "Install with: pip install scikit-learn"
            )

        self.n_folds = int(n_folds)
        self.n_estimators = int(n_estimators)

        self.mean_importance_: Optional[np.ndarray] = None
        self.std_importance_: Optional[np.ndarray] = None
        self.feature_names_: Optional[List[str]] = None

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        feature_names: Optional[List[str]] = None,
        random_state: Optional[int] = None,
    ) -> "CVFeatureImportance":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=int)

        if X.ndim != 2:
            raise ValueError("X must be 2D (n_samples, n_features)")
        if y.ndim != 1:
            raise ValueError("y must be 1D (n_samples,)")
        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y must have matching n_samples")

        n_samples, n_features = X.shape
        self.feature_names_ = (
            list(feature_names)
            if feature_names is not None
            else [f"feature_{i}" for i in range(n_features)]
        )

        kf = KFold(n_splits=min(self.n_folds, n_samples), shuffle=True, random_state=random_state)
        fold_importances: List[np.ndarray] = []

        for train_idx, _ in kf.split(X, y):
            model = RandomForestClassifier(
                n_estimators=self.n_estimators,
                random_state=random_state,
                n_jobs=-1,
            )
            model.fit(X[train_idx], y[train_idx])
            imp = np.asarray(getattr(model, "feature_importances_", np.zeros((n_features,))), dtype=float)
            if imp.shape[0] != n_features:
                continue
            fold_importances.append(imp)

        if not fold_importances:
            self.mean_importance_ = np.zeros((n_features,), dtype=float)
            self.std_importance_ = np.zeros((n_features,), dtype=float)
            return self

        mat = np.vstack(fold_importances)
        self.mean_importance_ = np.mean(mat, axis=0)
        self.std_importance_ = np.std(mat, axis=0)
        return self

    def get_ranked_features(self) -> List[Tuple[str, float, float]]:
        if self.mean_importance_ is None or self.std_importance_ is None or self.feature_names_ is None:
            raise ValueError("CVFeatureImportance must be fit() before get_ranked_features()")

        rows = [
            (self.feature_names_[i], float(self.mean_importance_[i]), float(self.std_importance_[i]))
            for i in range(len(self.feature_names_))
        ]
        rows.sort(key=lambda t: t[1], reverse=True)
        return rows


__all__ = ["CVFeatureImportance"]


