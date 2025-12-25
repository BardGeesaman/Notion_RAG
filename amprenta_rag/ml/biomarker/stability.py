"""Stability selection for biomarker discovery (bootstrap + LassoCV)."""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

try:
    from sklearn.linear_model import LassoCV
    from sklearn.preprocessing import StandardScaler

    SKLEARN_AVAILABLE = True
except Exception:  # noqa: BLE001
    SKLEARN_AVAILABLE = False


class StabilitySelector:
    """
    Stability selection using bootstrap resampling + L1 regularization (LassoCV).

    This implementation follows a pragmatic pattern used elsewhere in the repo:
    - bootstrap resample samples
    - standardize features
    - fit LassoCV
    - count non-zero coefficient selection frequency
    """

    def __init__(self, n_bootstrap: int = 50, threshold: float = 0.6):
        if not SKLEARN_AVAILABLE:
            raise ImportError(
                "scikit-learn is required for StabilitySelector (LassoCV). "
                "Install with: pip install scikit-learn"
            )

        self.n_bootstrap = int(n_bootstrap)
        self.threshold = float(threshold)

        self.selection_frequencies_: Optional[np.ndarray] = None
        self.mean_coefs_: Optional[np.ndarray] = None
        self.feature_names_: Optional[List[str]] = None

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        feature_names: Optional[List[str]] = None,
        random_state: Optional[int] = None,
    ) -> "StabilitySelector":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)

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

        selection_counts = np.zeros((n_features,), dtype=float)
        coef_sums = np.zeros((n_features,), dtype=float)

        rng = np.random.RandomState(random_state)

        for _ in range(self.n_bootstrap):
            idx = rng.choice(n_samples, size=n_samples, replace=True)
            Xb = X[idx, :]
            yb = y[idx]

            scaler = StandardScaler()
            Xb_scaled = scaler.fit_transform(Xb)

            model = LassoCV(cv=5, random_state=random_state, max_iter=10_000)
            try:
                model.fit(Xb_scaled, yb)
            except Exception:
                continue

            coefs = np.asarray(model.coef_, dtype=float)
            selected = np.abs(coefs) > 1e-8
            selection_counts[selected] += 1.0
            coef_sums += coefs

        self.selection_frequencies_ = selection_counts / float(max(self.n_bootstrap, 1))
        self.mean_coefs_ = coef_sums / float(max(self.n_bootstrap, 1))
        return self

    def get_ranked_features(self) -> List[Tuple[str, float, float]]:
        if self.selection_frequencies_ is None or self.mean_coefs_ is None or self.feature_names_ is None:
            raise ValueError("StabilitySelector must be fit() before get_ranked_features()")

        rows: List[Tuple[str, float, float]] = []
        for i, name in enumerate(self.feature_names_):
            freq = float(self.selection_frequencies_[i])
            if freq < self.threshold:
                continue
            rows.append((name, freq, float(self.mean_coefs_[i])))

        rows.sort(key=lambda t: t[1], reverse=True)
        return rows

