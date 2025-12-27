"""Tests for HTS plate full API endpoint."""

from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_plate_full_endpoint(monkeypatch) -> None:
    """Test GET /api/v1/hts/campaigns/{id}/plate/full endpoint."""
    from amprenta_rag.api.routers import hts as hts_router
    
    # Create test campaign ID
    test_campaign_id = uuid4()
    
    # Mock campaign exists check
    def mock_campaign_exists(campaign_id):
        if campaign_id != test_campaign_id:
            from fastapi import HTTPException
            raise HTTPException(status_code=404, detail="HTS campaign not found")
    
    # Mock get_plate_heatmap_data function
    def mock_get_plate_heatmap_data(campaign_id, include_all=False):
        # Verify that include_all=True is passed
        assert include_all is True, "include_all should be True for full plate endpoint"
        
        # Return mock well data
        return [
            SimpleNamespace(
                well_position="A01",
                normalized_value=0.8,
                z_score=2.5,
                hit_flag=True,
                compound_id=uuid4(),
                result_id="result1",
                asdict=lambda: {
                    "well_position": "A01",
                    "normalized_value": 0.8,
                    "z_score": 2.5,
                    "hit_flag": True,
                    "compound_id": test_campaign_id,  # Use test ID for consistency
                    "result_id": "result1"
                }
            ),
            SimpleNamespace(
                well_position="A02",
                normalized_value=0.2,
                z_score=0.5,
                hit_flag=False,
                compound_id=uuid4(),
                result_id="result2",
                asdict=lambda: {
                    "well_position": "A02",
                    "normalized_value": 0.2,
                    "z_score": 0.5,
                    "hit_flag": False,
                    "compound_id": test_campaign_id,
                    "result_id": "result2"
                }
            )
        ]
    
    # Apply mocks
    monkeypatch.setattr(hts_router, "_campaign_exists", mock_campaign_exists)
    monkeypatch.setattr(hts_router, "get_plate_heatmap_data", mock_get_plate_heatmap_data)
    
    # Test the endpoint
    resp = client.get(f"/api/v1/hts/campaigns/{test_campaign_id}/plate/full")
    assert resp.status_code == 200
    
    data = resp.json()
    assert isinstance(data, list)
    assert len(data) == 2
    
    # Verify response structure
    well1 = data[0]
    assert "well_position" in well1
    assert "normalized_value" in well1
    assert "z_score" in well1
    assert "hit_flag" in well1
    assert "compound_id" in well1
    assert "result_id" in well1
    
    # Verify data values
    assert well1["well_position"] == "A01"
    assert well1["normalized_value"] == 0.8
    assert well1["z_score"] == 2.5
    assert well1["hit_flag"] is True
    assert well1["result_id"] == "result1"
    
    well2 = data[1]
    assert well2["well_position"] == "A02"
    assert well2["normalized_value"] == 0.2
    assert well2["z_score"] == 0.5
    assert well2["hit_flag"] is False
    assert well2["result_id"] == "result2"


def test_plate_full_endpoint_not_found() -> None:
    """Test plate full endpoint with invalid campaign ID."""
    from amprenta_rag.api.routers import hts as hts_router
    from amprenta_rag.database.base import get_db
    
    # Mock database to return no campaign
    def mock_get_db():
        db = SimpleNamespace()
        
        def mock_query(model):
            query_obj = SimpleNamespace()
            
            def mock_filter(condition):
                filtered_obj = SimpleNamespace()
                filtered_obj.first = lambda: None  # No campaign found
                return filtered_obj
            
            query_obj.filter = mock_filter
            return query_obj
        
        db.query = mock_query
        yield db
    
    app.dependency_overrides[get_db] = mock_get_db
    
    try:
        invalid_campaign_id = uuid4()
        resp = client.get(f"/api/v1/hts/campaigns/{invalid_campaign_id}/plate/full")
        assert resp.status_code == 404
        assert "HTS campaign not found" in resp.json()["detail"]
    finally:
        # Clean up dependency override
        if get_db in app.dependency_overrides:
            del app.dependency_overrides[get_db]


def test_plate_full_endpoint_empty_data(monkeypatch) -> None:
    """Test plate full endpoint with campaign that has no well data."""
    from amprenta_rag.api.routers import hts as hts_router
    
    test_campaign_id = uuid4()
    
    # Mock campaign exists
    def mock_campaign_exists(campaign_id):
        pass  # Campaign exists, no exception
    
    # Mock get_plate_heatmap_data to return empty list
    def mock_get_plate_heatmap_data(campaign_id, include_all=False):
        assert include_all is True
        return []  # No well data
    
    monkeypatch.setattr(hts_router, "_campaign_exists", mock_campaign_exists)
    monkeypatch.setattr(hts_router, "get_plate_heatmap_data", mock_get_plate_heatmap_data)
    
    resp = client.get(f"/api/v1/hts/campaigns/{test_campaign_id}/plate/full")
    assert resp.status_code == 200
    
    data = resp.json()
    assert isinstance(data, list)
    assert len(data) == 0


def test_plate_full_endpoint_vs_regular_plate(monkeypatch) -> None:
    """Test that plate/full returns all wells while plate returns hits only."""
    from amprenta_rag.api.routers import hts as hts_router
    
    test_campaign_id = uuid4()
    
    # Mock campaign exists
    def mock_campaign_exists(campaign_id):
        pass
    
    # Track calls to get_plate_heatmap_data
    calls_made = []
    
    def mock_get_plate_heatmap_data(campaign_id, include_all=False):
        calls_made.append(include_all)
        
        if include_all:
            # Return all wells for /plate/full
            return [
                SimpleNamespace(asdict=lambda: {"well_position": "A01", "hit_flag": True}),
                SimpleNamespace(asdict=lambda: {"well_position": "A02", "hit_flag": False}),
            ]
        else:
            # Return hits only for /plate
            return [
                SimpleNamespace(asdict=lambda: {"well_position": "A01", "hit_flag": True}),
            ]
    
    monkeypatch.setattr(hts_router, "_campaign_exists", mock_campaign_exists)
    monkeypatch.setattr(hts_router, "get_plate_heatmap_data", mock_get_plate_heatmap_data)
    
    # Test regular plate endpoint (hits only)
    resp_plate = client.get(f"/api/v1/hts/campaigns/{test_campaign_id}/plate")
    assert resp_plate.status_code == 200
    plate_data = resp_plate.json()
    assert len(plate_data) == 1  # Only hits
    
    # Test full plate endpoint (all wells)
    resp_full = client.get(f"/api/v1/hts/campaigns/{test_campaign_id}/plate/full")
    assert resp_full.status_code == 200
    full_data = resp_full.json()
    assert len(full_data) == 2  # All wells
    
    # Verify correct parameters were passed
    assert len(calls_made) == 2
    assert calls_made[0] is False  # Regular plate endpoint
    assert calls_made[1] is True   # Full plate endpoint


def test_plate_full_endpoint_invalid_uuid() -> None:
    """Test plate full endpoint with malformed UUID."""
    resp = client.get("/api/v1/hts/campaigns/not-a-uuid/plate/full")
    assert resp.status_code == 422  # Validation error for invalid UUID format
