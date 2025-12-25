from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_dose_response_fit_endpoint(monkeypatch):
    from amprenta_rag.api.routers import explorer as explorer_router

    class FakeFit:
        model = "4PL"
        params = {"ec50": 10.0, "hill_slope": 1.0, "top": 1.0, "bottom": 0.1, "r_squared": 0.9}
        r_squared = 0.9

    monkeypatch.setattr(explorer_router, "fit_dose_response", lambda concentrations, responses, model="4PL": FakeFit())
    monkeypatch.setattr(explorer_router, "generate_fit_curve", lambda params, x_min, x_max, n_points=100: ([1, 10], [0.2, 0.8]))
    monkeypatch.setattr(explorer_router, "bootstrap_ci", lambda conc, resp, fit_func, n_bootstrap=500: {"ec50_ci": (5.0, 20.0)})

    resp = client.post(
        "/api/explorer/dose-response/fit",
        json={"concentrations": [1, 10, 100, 1000], "responses": [0.1, 0.2, 0.8, 0.95], "model": "4PL"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["ec50"] == 10.0
    assert data["curve_x"] == [1, 10]


def test_dose_response_compare_endpoint(monkeypatch):
    from amprenta_rag.api.routers import explorer as explorer_router

    class FakeFit:
        model = "4PL"
        params = {"ec50": 10.0, "hill_slope": 1.0, "top": 1.0, "bottom": 0.1}
        r_squared = 0.9

    monkeypatch.setattr(explorer_router, "fit_dose_response", lambda concentrations, responses, model="4PL": FakeFit())
    monkeypatch.setattr(explorer_router, "generate_fit_curve", lambda params, x_min, x_max, n_points=100: ([1, 10], [0.2, 0.8]))
    monkeypatch.setattr(explorer_router, "bootstrap_ci", lambda conc, resp, fit_func, n_bootstrap=500: {"ec50_ci": (5.0, 20.0)})

    payload = {
        "fits": [
            {"concentrations": [1, 10, 100, 1000], "responses": [0.1, 0.2, 0.8, 0.95], "model": "4PL"},
            {"concentrations": [1, 10, 100, 1000], "responses": [0.2, 0.3, 0.7, 0.9], "model": "4PL"},
        ]
    }
    resp = client.post("/api/explorer/dose-response/compare", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["results"]) == 2
    assert "comparison_table" in data


def test_timeseries_analyze_endpoint(monkeypatch):
    from amprenta_rag.api.routers import explorer as explorer_router

    monkeypatch.setattr(explorer_router, "analyze_timeseries", lambda values, timepoints: type("T", (), {"slope": 1.0, "pvalue": 0.01, "direction": "increasing"})())
    monkeypatch.setattr(explorer_router, "smooth_timeseries", lambda values, method="savgol", window=5: list(values))
    monkeypatch.setattr(explorer_router, "detect_changepoints", lambda values, threshold=2.0: [2])

    resp = client.post(
        "/api/explorer/timeseries/analyze",
        json={"values": [1, 2, 3], "timepoints": [0, 1, 2], "smooth_method": "savgol", "detect_changepoints": True},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["direction"] == "increasing"
    assert data["changepoint_indices"] == [2]


def test_timeseries_compare_endpoint(monkeypatch):
    from amprenta_rag.api.routers import explorer as explorer_router

    monkeypatch.setattr(explorer_router, "analyze_timeseries", lambda values, timepoints: type("T", (), {"slope": 1.0, "pvalue": 0.01, "direction": "increasing"})())
    monkeypatch.setattr(explorer_router, "compare_trajectories", lambda series_list: {"comparison_table": [{"label": "a"}], "cluster_labels": [0, 1]})

    resp = client.post(
        "/api/explorer/timeseries/compare",
        json={"series": [{"values": [1, 2, 3], "timepoints": [0, 1, 2]}, {"values": [3, 2, 1], "timepoints": [0, 1, 2]}]},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["results"]) == 2
    assert data["cluster_labels"] == [0, 1]


def test_invalid_data_returns_error():
    # Mismatched lengths should return warnings + NaN ec50, but still 200 (graceful handling).
    resp = client.post(
        "/api/explorer/dose-response/fit",
        json={"concentrations": [1, 10], "responses": [0.1], "model": "4PL"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["warnings"]
    assert float(data["ec50"]) == 0.0


