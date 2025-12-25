from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def _tanimoto(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Generalized Tanimoto similarity for non-negative vectors.

    For binary fingerprints `a` and fractional centroid `b`, this is a standard
    extension: sim = (a·b) / (||a||^2 + ||b||^2 - a·b).
    """

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if b.ndim != 1:
        raise ValueError("Centroid must be 1D")
    if a.ndim != 2:
        raise ValueError("X must be 2D")

    dot = a @ b
    a2 = np.sum(a * a, axis=1)
    b2 = float(np.sum(b * b))
    denom = a2 + b2 - dot
    denom = np.maximum(denom, 1e-12)
    sim = dot / denom
    return np.clip(sim, 0.0, 1.0)


class ApplicabilityChecker:
    def __init__(self, threshold: float = 0.3):
        self.threshold = float(threshold)
        self.training_centroid: Optional[np.ndarray] = None

    def fit(self, X_train: np.ndarray) -> "ApplicabilityChecker":
        X_train = np.asarray(X_train, dtype=float)
        if X_train.ndim != 2 or X_train.shape[0] < 1:
            raise ValueError("X_train must be 2D with at least 1 row")
        self.training_centroid = X_train.mean(axis=0).astype(np.float32, copy=False)
        return self

    def check(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if self.training_centroid is None:
            raise ValueError("ApplicabilityChecker not fit")
        X = np.asarray(X, dtype=float)
        sims = _tanimoto(X, self.training_centroid)
        in_domain = sims >= self.threshold
        return in_domain.astype(bool), sims.astype(float)

    def widen_ci(
        self, ci_low: np.ndarray, ci_high: np.ndarray, in_domain: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        ci_low = np.asarray(ci_low, dtype=float)
        ci_high = np.asarray(ci_high, dtype=float)
        in_domain = np.asarray(in_domain, dtype=bool)
        if ci_low.shape != ci_high.shape:
            raise ValueError("ci_low and ci_high must have same shape")
        if ci_low.shape[0] != in_domain.shape[0]:
            raise ValueError("in_domain length must match CI arrays")

        width = ci_high - ci_low
        # double width for OOD (i.e., not in_domain)
        widen = (~in_domain).astype(float)
        new_low = ci_low - 0.5 * width * widen
        new_high = ci_high + 0.5 * width * widen
        return new_low, new_high



