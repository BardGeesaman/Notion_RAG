from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.fixture()
def client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_get_targets_ok_calls_service(client, monkeypatch):
    from amprenta_rag.api.routers import sar as sar_router

    fn = Mock(return_value=[{"target": "CDK2", "compound_count": 2}])
    monkeypatch.setattr(sar_router.service, "list_targets", fn)

    resp = client.get("/api/v1/sar/targets?limit=10")
    assert resp.status_code == 200
    assert resp.json() == [{"target": "CDK2", "compound_count": 2}]
    fn.assert_called_once_with(limit=10)


@pytest.mark.api
def test_get_targets_limit_validation(client):
    resp = client.get("/api/v1/sar/targets?limit=0")
    assert resp.status_code == 422


@pytest.mark.api
def test_get_compounds_by_target_ok_calls_service(client, monkeypatch):
    from amprenta_rag.api.routers import sar as sar_router

    fn = Mock(
        return_value=[
            {
                "compound_id": "C1",
                "smiles": "CCO",
                "ic50": 10.0,
                "units": "nM",
                "assay_name": "A",
                "result_id": "R1",
            }
        ]
    )
    monkeypatch.setattr(sar_router.service, "get_compounds_by_target", fn)

    resp = client.get("/api/v1/sar/targets/CDK2/compounds?limit=123")
    assert resp.status_code == 200
    assert resp.json() == fn.return_value
    fn.assert_called_once_with(target="CDK2", limit=123)


@pytest.mark.api
def test_get_cliffs_ok_calls_service(client, monkeypatch):
    from amprenta_rag.api.routers import sar as sar_router

    fn = Mock(
        return_value=[
            {
                "compound_1": "C1",
                "smiles_1": "S1",
                "activity_1": 100.0,
                "compound_2": "C2",
                "smiles_2": "S2",
                "activity_2": 5.0,
                "similarity": 0.75,
                "fold_change": 20.0,
                "assay_id": None,
            }
        ]
    )
    monkeypatch.setattr(sar_router.service, "get_activity_cliffs_for_target", fn)

    resp = client.get(
        "/api/v1/sar/targets/CDK2/cliffs?similarity_threshold=0.7&fold_change=10&limit=5"
    )
    assert resp.status_code == 200
    assert resp.json() == fn.return_value
    fn.assert_called_once_with(
        target="CDK2",
        similarity_threshold=0.7,
        fold_change=10.0,
        limit=5,
    )


@pytest.mark.api
def test_get_cliffs_validation(client):
    resp = client.get("/api/v1/sar/targets/CDK2/cliffs?similarity_threshold=1.5")
    assert resp.status_code == 422
    resp = client.get("/api/v1/sar/targets/CDK2/cliffs?fold_change=0.5")
    assert resp.status_code == 422
    resp = client.get("/api/v1/sar/targets/CDK2/cliffs?limit=0")
    assert resp.status_code == 422


