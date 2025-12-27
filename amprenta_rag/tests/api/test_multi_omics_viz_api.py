"""Tests for multi-omics visualization API endpoints."""

from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_get_datasets_endpoint(monkeypatch):
    """Test GET /api/multi-omics-viz/datasets endpoint."""
    from amprenta_rag.api.routers import multi_omics_viz as mov_router
    
    # Mock database query
    mock_datasets = [
        SimpleNamespace(id=uuid4(), name="Dataset A", omics_type="transcriptomics"),
        SimpleNamespace(id=uuid4(), name="Dataset B", omics_type="proteomics"),
        SimpleNamespace(id=uuid4(), name="Dataset C", omics_type="metabolomics"),
    ]
    
    def mock_get_db():
        mock_db = SimpleNamespace()
        mock_db.query = lambda *args: SimpleNamespace(all=lambda: mock_datasets)
        yield mock_db
    
    app.dependency_overrides[mov_router.get_db] = mock_get_db
    
    resp = client.get("/api/multi-omics-viz/datasets")
    assert resp.status_code == 200
    data = resp.json()
    
    assert len(data) == 3
    assert all("id" in item for item in data)
    assert all("name" in item for item in data)
    assert all("omics_type" in item for item in data)
    assert data[0]["omics_type"] == "transcriptomics"
    assert data[1]["omics_type"] == "proteomics"
    assert data[2]["omics_type"] == "metabolomics"
    
    # Clean up dependency override
    del app.dependency_overrides[mov_router.get_db]


def test_alluvial_endpoint(monkeypatch):
    """Test POST /api/multi-omics-viz/alluvial endpoint."""
    from amprenta_rag.api.routers import multi_omics_viz as mov_router
    
    # Mock compute_alluvial_data
    mock_result = {
        "nodes": [
            {"id": "dataset_1", "label": "Dataset A", "color": "#E74C3C"},
            {"id": "transcriptomics", "label": "transcriptomics", "color": "#E74C3C"},
            {"id": "feature_tp53", "label": "TP53", "color": "#9E9E9E"},
        ],
        "links": [
            {"source": 0, "target": 1, "value": 100},
            {"source": 1, "target": 2, "value": 50},
        ]
    }
    
    def mock_compute_alluvial(dataset_ids, db):
        return mock_result
    
    monkeypatch.setattr(mov_router, "compute_alluvial_data", mock_compute_alluvial)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[mov_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/multi-omics-viz/alluvial",
            json={"dataset_ids": [str(uuid4()), str(uuid4())]}
        )
        assert resp.status_code == 200
        data = resp.json()
        
        # Verify AlluvialResponse schema
        assert "nodes" in data
        assert "links" in data
        assert len(data["nodes"]) == 3
        assert len(data["links"]) == 2
        
        # Verify node structure
        node = data["nodes"][0]
        assert "id" in node
        assert "label" in node
        assert "color" in node
        
        # Verify link structure
        link = data["links"][0]
        assert "source" in link
        assert "target" in link
        assert "value" in link
        
    finally:
        # Clean up dependency override
        del app.dependency_overrides[mov_router.get_db]


def test_upset_endpoint(monkeypatch):
    """Test POST /api/multi-omics-viz/upset endpoint."""
    from amprenta_rag.api.routers import multi_omics_viz as mov_router
    
    # Mock compute_upset_data
    mock_result = {
        "sets": [
            {"id": str(uuid4()), "name": "Dataset A", "omics_type": "transcriptomics", "color": "#E74C3C", "size": 100},
            {"id": str(uuid4()), "name": "Dataset B", "omics_type": "proteomics", "color": "#3498DB", "size": 80},
        ],
        "intersections": [
            {"sets": ["dataset_1", "dataset_2"], "count": 25, "bitmask": 3},
            {"sets": ["dataset_1"], "count": 75, "bitmask": 1},
            {"sets": ["dataset_2"], "count": 55, "bitmask": 2},
        ],
        "matrix": [
            {"key": "tp53", "presence": [1, 1], "label": "TP53", "dataset_count": 2},
            {"key": "brca1", "presence": [1, 0], "label": "BRCA1", "dataset_count": 1},
            {"key": "akt1", "presence": [0, 1], "label": "AKT1", "dataset_count": 1},
        ]
    }
    
    def mock_compute_upset(dataset_ids, db):
        return mock_result
    
    monkeypatch.setattr(mov_router, "compute_upset_data", mock_compute_upset)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[mov_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/multi-omics-viz/upset",
            json={"dataset_ids": [str(uuid4()), str(uuid4())]}
        )
        assert resp.status_code == 200
        data = resp.json()
        
        # Verify UpSetResponse schema
        assert "sets" in data
        assert "intersections" in data
        assert "matrix" in data
        
        # Verify sets structure
        assert len(data["sets"]) == 2
        set_item = data["sets"][0]
        assert "id" in set_item
        assert "name" in set_item
        assert "omics_type" in set_item
        assert "color" in set_item
        assert "size" in set_item
        
        # Verify intersections structure
        assert len(data["intersections"]) == 3
        intersection = data["intersections"][0]
        assert "sets" in intersection
        assert "count" in intersection
        assert "bitmask" in intersection
        
        # Verify matrix structure
        assert len(data["matrix"]) == 3
        assert data["matrix"][0]["key"] == "tp53"
        assert data["matrix"][0]["presence"] == [1, 1]
        
    finally:
        # Clean up dependency override
        del app.dependency_overrides[mov_router.get_db]


def test_alluvial_endpoint_invalid_uuid():
    """Test alluvial endpoint with invalid UUID."""
    resp = client.post(
        "/api/multi-omics-viz/alluvial",
        json={"dataset_ids": ["invalid-uuid"]}
    )
    assert resp.status_code == 400
    assert "Invalid UUID" in resp.json()["detail"]


def test_upset_endpoint_invalid_uuid():
    """Test upset endpoint with invalid UUID."""
    resp = client.post(
        "/api/multi-omics-viz/upset",
        json={"dataset_ids": ["not-a-uuid", "also-not-a-uuid"]}
    )
    assert resp.status_code == 400
    assert "Invalid UUID" in resp.json()["detail"]


def test_alluvial_endpoint_empty_dataset_ids(monkeypatch):
    """Test alluvial endpoint with empty dataset_ids."""
    from amprenta_rag.api.routers import multi_omics_viz as mov_router
    
    # Mock compute_alluvial_data to return empty result
    def mock_compute_alluvial(dataset_ids, db):
        return {"nodes": [], "links": []}
    
    monkeypatch.setattr(mov_router, "compute_alluvial_data", mock_compute_alluvial)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[mov_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/multi-omics-viz/alluvial",
            json={"dataset_ids": []}
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data == {"nodes": [], "links": []}
        
    finally:
        # Clean up dependency override
        del app.dependency_overrides[mov_router.get_db]
