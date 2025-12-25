from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np


def _require_sklearn():
    try:
        from sklearn.isotonic import IsotonicRegression  # type: ignore
        from sklearn.linear_model import LogisticRegression  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise ImportError("scikit-learn is required for calibration (pip install scikit-learn)") from e
    return IsotonicRegression, LogisticRegression


@dataclass(frozen=True)
class ReliabilityDiagram:
    bin_edges: List[float]
    bin_accs: List[float]
    bin_confs: List[float]
    ece: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "bin_edges": self.bin_edges,
            "bin_accs": self.bin_accs,
            "bin_confs": self.bin_confs,
            "ece": float(self.ece),
        }


class CalibrationWrapper:
    def __init__(self, method: str = "isotonic"):
        if method not in ("isotonic", "platt"):
            raise ValueError("method must be 'isotonic' or 'platt'")
        self.method = method
        self.calibrator: Optional[Any] = None

    def fit(self, y_pred: np.ndarray, y_true: np.ndarray) -> "CalibrationWrapper":
        y_pred = np.asarray(y_pred, dtype=float).reshape(-1)
        y_true = np.asarray(y_true, dtype=float).reshape(-1)
        if len(y_pred) != len(y_true):
            raise ValueError("y_pred and y_true must have same length")
        if len(y_pred) < 10:
            raise ValueError("Need at least 10 samples to fit calibrator")

        IsotonicRegression, LogisticRegression = _require_sklearn()

        if self.method == "isotonic":
            iso = IsotonicRegression(out_of_bounds="clip")
            iso.fit(y_pred, y_true)
            self.calibrator = iso
        else:
            # Platt scaling: logistic regression on single logit-like feature
            lr = LogisticRegression(solver="lbfgs")
            lr.fit(y_pred.reshape(-1, 1), y_true.astype(int))
            self.calibrator = lr
        return self

    def calibrate(self, y_pred: np.ndarray) -> np.ndarray:
        if self.calibrator is None:
            raise ValueError("Calibrator not fit")
        y_pred = np.asarray(y_pred, dtype=float).reshape(-1)

        if self.method == "isotonic":
            out = self.calibrator.transform(y_pred)
        else:
            out = self.calibrator.predict_proba(y_pred.reshape(-1, 1))[:, 1]
        return np.clip(np.asarray(out, dtype=float), 0.0, 1.0)

    def compute_ece(self, y_pred: np.ndarray, y_true: np.ndarray, n_bins: int = 10) -> float:
        rd = reliability_diagram(y_pred, y_true, n_bins=n_bins)
        return float(rd["ece"])


def reliability_diagram(y_pred: np.ndarray, y_true: np.ndarray, n_bins: int = 10) -> dict:
    y_pred = np.asarray(y_pred, dtype=float).reshape(-1)
    y_true = np.asarray(y_true, dtype=float).reshape(-1)
    if len(y_pred) != len(y_true):
        raise ValueError("y_pred and y_true must have same length")
    if n_bins <= 1:
        raise ValueError("n_bins must be > 1")

    y_pred = np.clip(y_pred, 0.0, 1.0)

    bin_edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    bin_accs: List[float] = []
    bin_confs: List[float] = []
    ece = 0.0

    for i in range(len(bin_edges) - 1):
        lo = bin_edges[i]
        hi = bin_edges[i + 1]
        if i == len(bin_edges) - 2:
            mask = (y_pred >= lo) & (y_pred <= hi)
        else:
            mask = (y_pred >= lo) & (y_pred < hi)

        if not mask.any():
            bin_accs.append(float("nan"))
            bin_confs.append(float("nan"))
            continue

        conf = float(np.mean(y_pred[mask]))
        acc = float(np.mean(y_true[mask]))
        bin_confs.append(conf)
        bin_accs.append(acc)
        ece += (mask.mean()) * abs(acc - conf)

    rd = ReliabilityDiagram(
        bin_edges=[float(x) for x in bin_edges.tolist()],
        bin_accs=[float(x) for x in bin_accs],
        bin_confs=[float(x) for x in bin_confs],
        ece=float(ece),
    )
    return rd.to_dict()



