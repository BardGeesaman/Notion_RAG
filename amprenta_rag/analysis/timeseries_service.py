"""Time-series analysis helpers (trend, smoothing, changepoints, forecasting).

This module wraps existing `amprenta_rag.analysis.timeseries.detect_trend` and adds
lightweight utilities commonly used in dashboards.

Examples:
    >>> analyze_timeseries([1, 2, 3], [0, 1, 2]).direction in {"increasing", "stable"}
    True
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Literal, Tuple

import numpy as np
from scipy.signal import savgol_filter
from scipy.stats import linregress

from amprenta_rag.analysis import timeseries


@dataclass(frozen=True)
class TimeseriesResult:
    slope: float
    pvalue: float
    direction: str


def analyze_timeseries(values: Iterable[float], timepoints: Iterable[float]) -> TimeseriesResult:
    """
    Analyze a single trajectory using the existing trend detector.

    Returns:
        TimeseriesResult(slope, pvalue, direction)
    """
    slope, p, direction = timeseries.detect_trend(values, timepoints)
    return TimeseriesResult(slope=float(slope), pvalue=float(p), direction=str(direction))


def smooth_timeseries(values: Iterable[float], method: Literal["savgol", "lowess"] = "savgol", window: int = 5) -> List[float]:
    """
    Smooth a time series.

    Methods:
      - savgol: Savitzky-Golay filter (SciPy)
      - lowess: LOWESS smoothing (statsmodels, optional)
    """
    y = np.asarray(list(values), dtype=float)
    if y.size == 0:
        return []
    m = (method or "savgol").lower()
    w = int(window)

    if m == "savgol":
        if w < 3:
            return list(map(float, y.tolist()))
        if w % 2 == 0:
            raise ValueError("window must be odd for savgol")
        w = min(w, y.size if y.size % 2 == 1 else y.size - 1)
        if w < 3:
            return list(map(float, y.tolist()))
        out = savgol_filter(y, window_length=w, polyorder=2, mode="interp")
        return list(map(float, out.tolist()))

    if m == "lowess":
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess  # type: ignore
        except Exception as e:  # noqa: BLE001
            raise ImportError("statsmodels is required for LOWESS smoothing (pip install statsmodels)") from e
        x = np.arange(y.size, dtype=float)
        frac = min(1.0, max(0.05, float(w) / float(y.size)))
        sm = lowess(y, x, frac=frac, return_sorted=False)
        return list(map(float, sm.tolist()))

    raise ValueError("method must be one of: savgol, lowess")


def detect_changepoints(values: Iterable[float], threshold: float = 2.0) -> List[int]:
    """
    Detect changepoints via rolling standard deviation spikes (heuristic).

    Returns indices (0-based) where rolling std exceeds `threshold`.
    """
    y = np.asarray(list(values), dtype=float)
    if y.size < 5:
        return []
    w = 5
    roll_std = np.array([np.nanstd(y[max(0, i - w + 1) : i + 1]) for i in range(y.size)], dtype=float)
    idx = np.where(roll_std > float(threshold))[0]
    return [int(i) for i in idx.tolist()]


def forecast_trend(values: Iterable[float], timepoints: Iterable[float], horizon: int = 5) -> Tuple[List[float], List[float]]:
    """
    Forecast future values via linear regression extrapolation.

    Returns (forecast_timepoints, forecast_values).
    """
    y = np.asarray(list(values), dtype=float)
    t = np.asarray(list(timepoints), dtype=float)
    if y.size != t.size or y.size < 2:
        return [], []

    mask = ~(np.isnan(y) | np.isnan(t))
    y = y[mask]
    t = t[mask]
    if y.size < 2:
        return [], []

    slope, intercept, *_ = linregress(t, y)
    h = int(horizon)
    if h < 1:
        return [], []

    # Step based on median delta.
    dt = np.diff(np.sort(t))
    step = float(np.nanmedian(dt)) if dt.size else 1.0
    if not np.isfinite(step) or step == 0:
        step = 1.0

    t_last = float(np.nanmax(t))
    ft = [t_last + step * (i + 1) for i in range(h)]
    fv = [float(intercept + slope * tt) for tt in ft]
    return ft, fv


def compare_trajectories(series_list: List[Any]) -> Dict[str, Any]:
    """
    Compare slopes across multiple series.

    Accepts a list of:
      - list/tuple/np.ndarray of values (timepoints inferred as 0..n-1), or
      - dict {"values": [...], "timepoints": [...], "label": "..."}.

    Returns:
      {"comparison_table": [...], "cluster_labels": [...]}
    """
    rows: List[Dict[str, Any]] = []
    slopes: List[float] = []

    for i, item in enumerate(series_list or []):
        label = f"series_{i}"
        if isinstance(item, dict):
            label = str(item.get("label") or label)
            vals = item.get("values") or []
            tps = item.get("timepoints")
        else:
            vals = item
            tps = None

        y = np.asarray(list(vals), dtype=float)
        if tps is None:
            t = np.arange(y.size, dtype=float)
        else:
            t = np.asarray(list(tps), dtype=float)

        s, p, d = timeseries.detect_trend(y.tolist(), t.tolist())
        slopes.append(float(s))
        rows.append({"label": label, "slope": float(s), "pvalue": float(p), "direction": str(d)})

    # Cluster labels (best-effort): KMeans on slopes if sklearn available; else sign buckets.
    cluster_labels: List[int] = []
    try:
        from sklearn.cluster import KMeans  # type: ignore

        arr = np.asarray(slopes, dtype=float).reshape(-1, 1)
        k = min(3, max(1, arr.shape[0]))
        if arr.shape[0] < 2:
            cluster_labels = [0 for _ in slopes]
        else:
            km = KMeans(n_clusters=k, n_init=10, random_state=42)
            cluster_labels = [int(x) for x in km.fit_predict(arr).tolist()]
    except Exception:
        for s in slopes:
            if not np.isfinite(s):
                cluster_labels.append(0)
            elif s > 0:
                cluster_labels.append(1)
            elif s < 0:
                cluster_labels.append(2)
            else:
                cluster_labels.append(0)

    return {"comparison_table": rows, "cluster_labels": cluster_labels}


__all__ = [
    "TimeseriesResult",
    "analyze_timeseries",
    "smooth_timeseries",
    "detect_changepoints",
    "forecast_trend",
    "compare_trajectories",
]


