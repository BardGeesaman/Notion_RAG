from __future__ import annotations

import pytest
from fastapi.testclient import TestClient
from uuid import uuid4

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_get_methods_returns_list():
    resp = client.get("/api/biomarker/methods")
    assert resp.status_code == 200
    data = resp.json()
    assert isinstance(data, list)
    assert set(data) >= {"statistical", "stability", "importance"}


def test_discover_endpoint_returns_consensus(monkeypatch):
    from amprenta_rag.api.routers import biomarker as biomarker_router

    class FakeSvc:
        def discover(self, experiment_id, group1, group2, methods=None):  # noqa: ANN001
            assert group1 == ["S1"]
            assert group2 == ["S2"]
            return {
                "consensus_ranking": [
                    {"feature": "A", "avg_rank": 1.0, "methods": 1},
                    {"feature": "B", "avg_rank": 2.0, "methods": 1},
                ],
                "method_results": {
                    "statistical": [
                        {"feature": "A", "t_stat": 1.0, "p_value": 0.01, "p_adj": 0.01},
                        {"feature": "B", "t_stat": 0.5, "p_value": 0.5, "p_adj": 0.5},
                    ]
                },
            }

    monkeypatch.setattr(biomarker_router, "BiomarkerDiscoveryService", lambda: FakeSvc())

    exp_id = uuid4()
    resp = client.post(
        "/api/biomarker/discover",
        json={
            "experiment_id": str(exp_id),
            "group1_samples": ["S1"],
            "group2_samples": ["S2"],
            "methods": ["statistical"],
            "fdr_threshold": 0.05,
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert "consensus_ranking" in data
    assert "method_results" in data
    # FDR filtering should keep only feature A
    assert [r["feature"] for r in data["consensus_ranking"]] == ["A"]
    assert [r["feature"] for r in data["method_results"]["statistical"]] == ["A"]


