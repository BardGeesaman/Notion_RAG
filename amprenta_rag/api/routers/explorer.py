"""Dose-response + time-series explorer API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
from fastapi import APIRouter

from amprenta_rag.analysis.dose_response_service import (
    bootstrap_ci,
    compare_curves,
    fit_3pl,
    fit_4pl,
    fit_dose_response,
    generate_fit_curve,
)
from amprenta_rag.analysis.timeseries_service import (
    analyze_timeseries,
    compare_trajectories,
    detect_changepoints,
    smooth_timeseries,
)
from amprenta_rag.api.schemas import (
    DoseResponseCompareRequest,
    DoseResponseCompareResponse,
    DoseResponseFitRequest,
    DoseResponseFitResponse,
    TimeseriesAnalyzeRequest,
    TimeseriesAnalyzeResponse,
    TrajectoryCompareRequest,
    TrajectoryCompareResponse,
)


router = APIRouter(prefix="/explorer", tags=["explorer"])


def _nan_fit(warnings: List[str]) -> DoseResponseFitResponse:
    return DoseResponseFitResponse(
        ec50=0.0,
        ec50_ci=None,
        hill_slope=0.0,
        top=0.0,
        bottom=None,
        r_squared=None,
        curve_x=[],
        curve_y=[],
        ci_lower=None,
        ci_upper=None,
        warnings=warnings,
    )


def _fit_one(req: DoseResponseFitRequest) -> DoseResponseFitResponse:
    warnings: List[str] = []
    try:
        concentrations = list(map(float, req.concentrations))
        responses = list(map(float, req.responses))
    except Exception:
        return _nan_fit(["Invalid numeric inputs"])

    model = (req.model or "4PL").strip()
    if len(concentrations) != len(responses):
        return _nan_fit(["concentrations and responses must have same length"])

    # Default curve bounds from data.
    try:
        x_min = float(np.nanmin(np.asarray(concentrations, dtype=float)))
        x_max = float(np.nanmax(np.asarray(concentrations, dtype=float)))
    except Exception:
        x_min, x_max = 1.0, 1.0

    try:
        fit = fit_dose_response(concentrations, responses, model=model)  # type: ignore[arg-type]
    except ImportError as e:
        return _nan_fit([str(e)])
    except Exception as e:  # noqa: BLE001
        return _nan_fit([f"Fit failed: {e}"])

    params = fit.params or {}
    ec50 = float(params.get("ec50", float("nan")))
    hill = float(params.get("hill_slope", float("nan")))
    top = float(params.get("top", float("nan")))
    bottom = params.get("bottom")
    r2 = fit.r_squared

    # Smooth curve from fitted params (best-effort)
    curve_x: List[float]
    curve_y: List[float]
    try:
        curve_x, curve_y = generate_fit_curve(params, x_min=x_min, x_max=x_max, n_points=100)
    except Exception as e:  # noqa: BLE001
        warnings.append(f"Curve generation failed: {e}")
        curve_x, curve_y = [], []

    ec50_ci = None
    if (req.model or "").strip().lower() != "bayesian_4pl":
        # Bootstrap CI for speed (API cap)
        try:
            base_fit = fit_4pl if str(fit.model).upper() == "4PL" else fit_3pl
            cis = bootstrap_ci(concentrations, responses, base_fit, n_bootstrap=500)
            ec50_ci = tuple(cis.get("ec50_ci")) if cis.get("ec50_ci") else None  # type: ignore[arg-type]
        except Exception as e:  # noqa: BLE001
            warnings.append(f"Bootstrap CI failed: {e}")
    else:
        if "ec50_ci" in params:
            try:
                lo, hi = params["ec50_ci"]
                ec50_ci = (float(lo), float(hi))
            except Exception:
                ec50_ci = None
        if not np.isfinite(top):
            warnings.append("Bayesian fit does not return top/bottom in this MVP response.")

    return DoseResponseFitResponse(
        ec50=ec50,
        ec50_ci=ec50_ci,
        hill_slope=hill,
        top=top,
        bottom=float(bottom) if bottom is not None else None,
        r_squared=float(r2) if r2 is not None else None,
        curve_x=curve_x,
        curve_y=curve_y,
        ci_lower=None,
        ci_upper=None,
        warnings=warnings,
    )


@router.post("/dose-response/fit", response_model=DoseResponseFitResponse)
def dose_response_fit(payload: DoseResponseFitRequest) -> DoseResponseFitResponse:
    return _fit_one(payload)


@router.post("/dose-response/compare", response_model=DoseResponseCompareResponse)
def dose_response_compare(payload: DoseResponseCompareRequest) -> DoseResponseCompareResponse:
    fits = payload.fits or []
    if len(fits) > 50:
        fits = fits[:50]
    results = [_fit_one(f) for f in fits]

    comp_inputs: List[Dict[str, Any]] = []
    for i, r in enumerate(results):
        comp_inputs.append({"label": f"fit_{i}", "ec50": r.ec50, "hill_slope": r.hill_slope, "r_squared": r.r_squared})
    comparison = compare_curves(comp_inputs)
    return DoseResponseCompareResponse(results=results, comparison_table=comparison)


@router.post("/timeseries/analyze", response_model=TimeseriesAnalyzeResponse)
def timeseries_analyze(payload: TimeseriesAnalyzeRequest) -> TimeseriesAnalyzeResponse:
    try:
        res = analyze_timeseries(payload.values, payload.timepoints)
        smoothed = None
        if payload.smooth_method:
            smoothed = smooth_timeseries(payload.values, method=str(payload.smooth_method), window=5)
        cps = None
        if payload.detect_changepoints:
            cps = detect_changepoints(payload.values, threshold=float(payload.changepoint_threshold))
        return TimeseriesAnalyzeResponse(
            slope=float(res.slope),
            pvalue=float(res.pvalue),
            direction=str(res.direction),
            smoothed_values=smoothed,
            changepoint_indices=cps,
        )
    except ImportError:
        return TimeseriesAnalyzeResponse(
            slope=0.0,
            pvalue=1.0,
            direction="stable",
            smoothed_values=None,
            changepoint_indices=None,
        )
    except Exception:  # noqa: BLE001
        return TimeseriesAnalyzeResponse(
            slope=0.0,
            pvalue=1.0,
            direction="stable",
            smoothed_values=None,
            changepoint_indices=None,
        )


@router.post("/timeseries/compare", response_model=TrajectoryCompareResponse)
def timeseries_compare(payload: TrajectoryCompareRequest) -> TrajectoryCompareResponse:
    series = payload.series or []
    if len(series) > 50:
        series = series[:50]

    results: List[TimeseriesAnalyzeResponse] = []
    compare_in: List[Dict[str, Any]] = []

    for i, s in enumerate(series):
        r = timeseries_analyze(s)
        results.append(r)
        compare_in.append({"values": s.values, "timepoints": s.timepoints, "label": f"series_{i}"})

    comp = compare_trajectories(compare_in)
    return TrajectoryCompareResponse(
        results=results,
        comparison_table={"comparison_table": comp.get("comparison_table")},
        cluster_labels=comp.get("cluster_labels"),
    )


__all__ = ["router"]


