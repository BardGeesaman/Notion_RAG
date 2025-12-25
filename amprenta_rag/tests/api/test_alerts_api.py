from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
rdkit = pytest.importorskip("rdkit")
_ = rdkit

from fastapi.testclient import TestClient  # noqa: E402

from amprenta_rag.api.main import app  # noqa: E402


client = TestClient(app)


def test_check_endpoint_returns_alerts():
    # Rhodanine-like PAINS
    resp = client.post("/api/alerts/check", json={"smiles": "O=C1NC(=S)SC1=Cc1ccccc1", "filters": ["pains"]})
    assert resp.status_code == 200
    data = resp.json()
    assert data["traffic_light"] == "RED"
    assert data["alert_count"] >= 1
    assert any(a["alert_type"] == "PAINS" for a in data["alerts"])


def test_check_endpoint_clean_compound():
    # Ibuprofen should be clean under these filters
    resp = client.post("/api/alerts/check", json={"smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["traffic_light"] == "GREEN"
    assert data["alert_count"] == 0


def test_batch_endpoint_multiple():
    payload = {
        "smiles_list": ["CC(=O)Cl", "O=[N+]([O-])c1ccccc1", "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"],
        "filters": ["lilly", "brenk"],
    }
    resp = client.post("/api/alerts/batch", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_checked"] == 3
    assert len(data["results"]) == 3


def test_batch_endpoint_limit_enforced():
    payload = {"smiles_list": ["CCO"] * 101}
    resp = client.post("/api/alerts/batch", json=payload)
    assert resp.status_code == 400


def test_filters_endpoint_lists_all():
    resp = client.get("/api/alerts/filters")
    assert resp.status_code == 200
    data = resp.json()
    names = {d["name"] for d in data}
    assert names == {"pains", "brenk", "lilly"}


