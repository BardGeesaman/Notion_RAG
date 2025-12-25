from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_list_targets_endpoint(monkeypatch):
    from amprenta_rag.api.routers import qsar as qsar_router

    class FakePred:
        def list_available_targets(self):  # noqa: ANN001
            return [
                {"target": "EGFR", "model_name": "qsar_EGFR_ensemble", "version": "1.0.0", "metrics": {"auc": 0.9}},
            ]

    monkeypatch.setattr(qsar_router, "TargetQSARPredictor", lambda: FakePred())
    resp = client.get("/api/qsar/targets")
    assert resp.status_code == 200
    data = resp.json()
    assert data[0]["target"] == "EGFR"


def test_predict_endpoint_returns_results(monkeypatch):
    from amprenta_rag.api.routers import qsar as qsar_router

    class FakePred:
        def predict(self, smiles_list, targets):  # noqa: ANN001
            assert targets == ["EGFR"]
            return [
                {
                    "smiles": smiles_list[0],
                    "predictions": {
                        "EGFR": {
                            "probability": 0.8,
                            "std": 0.1,
                            "ci_low": 0.6,
                            "ci_high": 1.0,
                            "in_domain": True,
                            "similarity": 0.8,
                            "calibrated": True,
                            "active": True,
                        }
                    },
                    "error": None,
                }
            ]

    monkeypatch.setattr(qsar_router, "TargetQSARPredictor", lambda: FakePred())
    resp = client.post("/api/qsar/predict", json={"smiles_list": ["CCO"], "targets": ["EGFR"]})
    assert resp.status_code == 200
    data = resp.json()
    assert data["results"][0]["predictions"]["EGFR"]["active"] is True


def test_predict_endpoint_limit_enforced(monkeypatch):
    # No predictor should be constructed if limit is exceeded.
    resp = client.post("/api/qsar/predict", json={"smiles_list": ["C"] * 101, "targets": ["EGFR"]})
    assert resp.status_code == 400


