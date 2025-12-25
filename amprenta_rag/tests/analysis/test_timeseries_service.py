from __future__ import annotations

import pytest


def test_analyze_timeseries_returns_result():
    from amprenta_rag.analysis.timeseries_service import analyze_timeseries

    out = analyze_timeseries([1, 2, 3, 4], [0, 1, 2, 3])
    assert hasattr(out, "slope")
    assert hasattr(out, "pvalue")
    assert hasattr(out, "direction")


def test_smooth_timeseries_savgol():
    from amprenta_rag.analysis.timeseries_service import smooth_timeseries

    vals = [1, 2, 3, 4, 5, 6, 7]
    sm = smooth_timeseries(vals, method="savgol", window=5)
    assert len(sm) == len(vals)


def test_smooth_timeseries_validates_odd_window():
    from amprenta_rag.analysis.timeseries_service import smooth_timeseries

    with pytest.raises(ValueError):
        smooth_timeseries([1, 2, 3, 4, 5], method="savgol", window=4)


def test_detect_changepoints_finds_spikes():
    from amprenta_rag.analysis.timeseries_service import detect_changepoints

    # Construct a short high-variance region to trigger rolling-std changepoints.
    vals = [1, 1, 1, 1, 10, -10, 10, -10, 1, 1, 1]
    cps = detect_changepoints(vals, threshold=2.0)
    assert cps


def test_forecast_trend_extrapolates():
    from amprenta_rag.analysis.timeseries_service import forecast_trend

    t, y = [0, 1, 2, 3], [1, 2, 3, 4]
    ft, fv = forecast_trend(y, t, horizon=3)
    assert len(ft) == 3
    assert len(fv) == 3
    assert fv[-1] > fv[0]


def test_compare_trajectories_returns_labels():
    from amprenta_rag.analysis.timeseries_service import compare_trajectories

    out = compare_trajectories(
        [
            {"label": "a", "values": [1, 2, 3], "timepoints": [0, 1, 2]},
            {"label": "b", "values": [3, 2, 1], "timepoints": [0, 1, 2]},
        ]
    )
    assert "comparison_table" in out
    assert "cluster_labels" in out
    assert len(out["cluster_labels"]) == 2


