"""Dose-response curve fitting services (3PL/4PL + bootstrap CIs).

This module provides frequentist curve fitting via SciPy and optional Bayesian fitting
via PyMC (delegating to `amprenta_rag.analysis.bayesian_dose_response`).

Examples:
    >>> conc = [0.1, 1, 10, 100]
    >>> resp = [0.05, 0.2, 0.8, 0.95]
    >>> out = fit_4pl(conc, resp)
    >>> round(out["ec50"], 3) > 0
    True
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Literal, Tuple

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import bootstrap


def _as_xy(conc: List[float] | np.ndarray, resp: List[float] | np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(conc, dtype=float).reshape(-1)
    y = np.asarray(resp, dtype=float).reshape(-1)
    if x.shape != y.shape:
        raise ValueError("concentrations and responses must have the same length")
    if x.size < 4:
        raise ValueError("Need at least 4 points for dose-response fit")
    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]
    y = y[mask]
    if x.size < 4:
        raise ValueError("Not enough valid (non-NaN) points for fit")
    if np.any(x <= 0):
        raise ValueError("All concentrations must be > 0 for logistic models")
    return x, y


def _r_squared(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = float(np.nansum((y_true - y_pred) ** 2))
    ss_tot = float(np.nansum((y_true - np.nanmean(y_true)) ** 2))
    if ss_tot <= 0:
        return float("nan")
    return 1.0 - (ss_res / ss_tot)


def _logistic_3pl(x: np.ndarray, top: float, ec50: float, hill: float) -> np.ndarray:
    # y = top / (1 + (x/ec50)^(-hill))
    x = np.asarray(x, dtype=float)
    return float(top) / (1.0 + (x / float(ec50)) ** (-float(hill)))


def _logistic_4pl(x: np.ndarray, bottom: float, top: float, ec50: float, hill: float) -> np.ndarray:
    # y = bottom + (top-bottom) / (1 + (x/ec50)^(-hill))
    x = np.asarray(x, dtype=float)
    b = float(bottom)
    t = float(top)
    return b + (t - b) / (1.0 + (x / float(ec50)) ** (-float(hill)))


def fit_3pl(conc: List[float] | np.ndarray, resp: List[float] | np.ndarray) -> Dict[str, float]:
    """
    Fit a 3-parameter logistic (3PL) model:
        y = top / (1 + (x/ec50)^(-hill))
    """
    x, y = _as_xy(conc, resp)
    top0 = float(np.nanmax(y))
    ec50_0 = float(np.nanmedian(x))
    hill0 = 1.0

    bounds = (
        [min(0.0, float(np.nanmin(y)) - 10.0), 1e-12, -10.0],  # top, ec50, hill
        [float(np.nanmax(y)) + 10.0, float(np.nanmax(x) * 1e6), 10.0],
    )
    popt, _pcov = curve_fit(
        lambda xx, top, ec50, hill: _logistic_3pl(xx, top, ec50, hill),
        x,
        y,
        p0=[top0, ec50_0, hill0],
        bounds=bounds,
        maxfev=20000,
    )
    top, ec50, hill = map(float, popt)
    yhat = _logistic_3pl(x, top=top, ec50=ec50, hill=hill)
    return {"ec50": ec50, "hill_slope": hill, "top": top, "r_squared": float(_r_squared(y, yhat))}


def fit_4pl(conc: List[float] | np.ndarray, resp: List[float] | np.ndarray) -> Dict[str, float]:
    """
    Fit a 4-parameter logistic (4PL) model:
        y = bottom + (top-bottom) / (1 + (x/ec50)^(-hill))
    """
    x, y = _as_xy(conc, resp)
    bottom0 = float(np.nanmin(y))
    top0 = float(np.nanmax(y))
    ec50_0 = float(np.nanmedian(x))
    hill0 = 1.0

    y_min = float(np.nanmin(y))
    y_max = float(np.nanmax(y))
    y_rng = max(1e-8, y_max - y_min)

    bounds = (
        [y_min - 10.0 * y_rng, y_min - 10.0 * y_rng, 1e-12, -10.0],
        [y_max + 10.0 * y_rng, y_max + 10.0 * y_rng, float(np.nanmax(x) * 1e6), 10.0],
    )
    popt, _pcov = curve_fit(
        lambda xx, bottom, top, ec50, hill: _logistic_4pl(xx, bottom, top, ec50, hill),
        x,
        y,
        p0=[bottom0, top0, ec50_0, hill0],
        bounds=bounds,
        maxfev=40000,
    )
    bottom, top, ec50, hill = map(float, popt)
    yhat = _logistic_4pl(x, bottom=bottom, top=top, ec50=ec50, hill=hill)
    return {
        "ec50": ec50,
        "hill_slope": hill,
        "top": top,
        "bottom": bottom,
        "r_squared": float(_r_squared(y, yhat)),
    }


@dataclass(frozen=True)
class DoseResponseFit:
    model: str
    params: Dict[str, Any]
    r_squared: float | None = None
    method: str = "scipy"


def fit_dose_response(
    concentrations: List[float],
    responses: List[float],
    model: Literal["3PL", "4PL", "bayesian_4pl"] = "4PL",
) -> DoseResponseFit:
    """
    Dispatcher for dose-response fitting.

    Models:
    - "3PL": SciPy curve_fit 3PL
    - "4PL": SciPy curve_fit 4PL
    - "bayesian_4pl": PyMC 4PL (optional dependency; wraps fit_bayesian_dose_response)
    """
    m = (model or "4PL").upper()
    if m == "3PL":
        out = fit_3pl(concentrations, responses)
        return DoseResponseFit(model="3PL", params=out, r_squared=out.get("r_squared"), method="scipy")
    if m == "4PL":
        out = fit_4pl(concentrations, responses)
        return DoseResponseFit(model="4PL", params=out, r_squared=out.get("r_squared"), method="scipy")
    if m in ("BAYESIAN_4PL", "BAYESIAN_4PL".upper()):
        from amprenta_rag.analysis.bayesian_dose_response import fit_bayesian_dose_response

        bout = fit_bayesian_dose_response(concentrations=list(map(float, concentrations)), responses=list(map(float, responses)))
        # Normalize field names into a params dict.
        params = {
            "ec50": float(bout.get("ec50_mean")),
            "ec50_ci": bout.get("ec50_ci"),
            "hill_slope": float(bout.get("hill_slope")),
        }
        return DoseResponseFit(model="bayesian_4pl", params=params, r_squared=None, method="pymc")
    raise ValueError("model must be one of: 3PL, 4PL, bayesian_4pl")


def generate_fit_curve(params: Dict[str, Any], x_min: float, x_max: float, n_points: int = 100) -> Tuple[List[float], List[float]]:
    """
    Generate a smooth fitted curve from fit parameters.

    Uses log-spaced x when x_min/x_max are positive (typical for concentrations).
    """
    if n_points < 2:
        raise ValueError("n_points must be >= 2")
    xmin = float(x_min)
    xmax = float(x_max)
    if xmin <= 0 or xmax <= 0:
        xs = np.linspace(xmin, xmax, int(n_points))
    else:
        xs = np.logspace(np.log10(xmin), np.log10(xmax), int(n_points))

    if "bottom" in params:
        ys = _logistic_4pl(xs, bottom=float(params["bottom"]), top=float(params["top"]), ec50=float(params["ec50"]), hill=float(params["hill_slope"]))
    else:
        ys = _logistic_3pl(xs, top=float(params["top"]), ec50=float(params["ec50"]), hill=float(params["hill_slope"]))
    return list(map(float, xs.tolist())), list(map(float, ys.tolist()))


def bootstrap_ci(
    conc: List[float],
    resp: List[float],
    fit_func: Callable[[List[float] | np.ndarray, List[float] | np.ndarray], Dict[str, float]],
    n_bootstrap: int = 1000,
) -> Dict[str, Tuple[float, float]]:
    """
    Estimate confidence intervals for EC50, Hill slope, and R^2 via bootstrap.

    Returns:
        {"ec50_ci": (low, high), "hill_ci": (low, high), "r_squared_ci": (low, high)}
    """
    x, y = _as_xy(conc, resp)

    def stat(xb, yb):  # noqa: ANN001
        # bootstrap passes resampled arrays; flatten for fit.
        xx = np.asarray(xb, dtype=float).reshape(-1)
        yy = np.asarray(yb, dtype=float).reshape(-1)
        out = fit_func(xx, yy)
        return np.array([out.get("ec50"), out.get("hill_slope"), out.get("r_squared")], dtype=float)

    res = bootstrap(
        data=(x, y),
        statistic=stat,
        paired=True,
        n_resamples=int(n_bootstrap),
        confidence_level=0.95,
        vectorized=False,
        random_state=42,
        method="percentile",
    )
    lo = res.confidence_interval.low
    hi = res.confidence_interval.high
    return {
        "ec50_ci": (float(lo[0]), float(hi[0])),
        "hill_ci": (float(lo[1]), float(hi[1])),
        "r_squared_ci": (float(lo[2]), float(hi[2])),
    }


def compare_curves(fits: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Compare EC50 values across multiple fit dicts.

    Input fits should include an 'ec50' key (and optionally 'label').
    """
    rows = []
    ec50s = []
    for i, f in enumerate(fits or []):
        if not isinstance(f, dict):
            continue
        label = f.get("label") or f.get("name") or f"fit_{i}"
        ec50 = f.get("ec50")
        try:
            ec50_v = float(ec50)
        except Exception:
            ec50_v = float("nan")
        ec50s.append(ec50_v)
        rows.append({"label": str(label), "ec50": ec50_v, "hill_slope": f.get("hill_slope"), "r_squared": f.get("r_squared")})

    base = ec50s[0] if ec50s else float("nan")
    fold = []
    for v in ec50s:
        if base and np.isfinite(base) and np.isfinite(v) and base > 0:
            fold.append(float(v / base))
        else:
            fold.append(float("nan"))
    return {"parameter_table": rows, "ec50_fold_changes": fold}


__all__ = [
    "DoseResponseFit",
    "fit_3pl",
    "fit_4pl",
    "fit_dose_response",
    "generate_fit_curve",
    "bootstrap_ci",
    "compare_curves",
]


