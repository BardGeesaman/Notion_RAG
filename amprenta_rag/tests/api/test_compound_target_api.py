from __future__ import annotations

import pytest
from fastapi.testclient import TestClient
from uuid import uuid4

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_network_endpoint_returns_nodes_edges(monkeypatch):
    from amprenta_rag.api.routers import compound_target as ct_router

    class FakeSvc:
        def get_compound_target_network(self, compound_ids=None, target_ids=None, filters=None, **kwargs):  # noqa: ANN001
            return {
                "nodes": [{"data": {"id": "compound:1"}}, {"data": {"id": "feature:2"}}],
                "edges": [{"data": {"id": "e1", "source": "compound:1", "target": "feature:2"}}],
                "meta": {"node_count": 2, "edge_count": 1},
            }

    monkeypatch.setattr(ct_router, "CompoundTargetNetworkService", lambda: FakeSvc())
    resp = client.post("/api/network/compound-target", json={"compound_ids": [], "target_ids": [], "filters": {}})
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["nodes"]) == 2
    assert len(data["edges"]) == 1


def test_network_endpoint_with_filters(monkeypatch):
    from amprenta_rag.api.routers import compound_target as ct_router

    seen = {}

    class FakeSvc:
        def get_compound_target_network(self, compound_ids=None, target_ids=None, filters=None, **kwargs):  # noqa: ANN001
            seen["filters"] = filters
            return {"nodes": [], "edges": [], "meta": {"node_count": 0, "edge_count": 0}}

    monkeypatch.setattr(ct_router, "CompoundTargetNetworkService", lambda: FakeSvc())
    resp = client.post(
        "/api/network/compound-target",
        json={
            "compound_ids": [],
            "target_ids": [],
            "filters": {"ic50_range": {"min_nm": 1.0, "max_nm": 100.0}, "activity_type": "IC50"},
        },
    )
    assert resp.status_code == 200
    assert seen["filters"]["activity_type"] == "IC50"


def test_expand_compound_endpoint(monkeypatch):
    from amprenta_rag.api.routers import compound_target as ct_router

    class FakeSvc:
        def expand_from_compound(self, compound_id):  # noqa: ANN001
            return [uuid4(), uuid4()]

    monkeypatch.setattr(ct_router, "CompoundTargetNetworkService", lambda: FakeSvc())
    cid = uuid4()
    resp = client.get(f"/api/network/compound-target/expand/compound/{cid}")
    assert resp.status_code == 200
    data = resp.json()
    assert "target_ids" in data
    assert len(data["target_ids"]) == 2


def test_expand_target_endpoint(monkeypatch):
    from amprenta_rag.api.routers import compound_target as ct_router

    class FakeSvc:
        def expand_from_target(self, target_id):  # noqa: ANN001
            return [uuid4()]

    monkeypatch.setattr(ct_router, "CompoundTargetNetworkService", lambda: FakeSvc())
    tid = uuid4()
    resp = client.get(f"/api/network/compound-target/expand/target/{tid}")
    assert resp.status_code == 200
    data = resp.json()
    assert data["compound_ids"]


def test_invalid_uuid_returns_422():
    resp = client.get("/api/network/compound-target/expand/compound/not-a-uuid")
    assert resp.status_code == 422


