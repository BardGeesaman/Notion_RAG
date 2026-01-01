"""API tests for lifecycle management endpoints."""

import pytest
from unittest.mock import patch, MagicMock
from uuid import uuid4

from fastapi.testclient import TestClient
from amprenta_rag.api.main import app


@pytest.fixture
def client():
    """Create test client with auth override."""
    mock_user = MagicMock()
    mock_user.id = uuid4()
    mock_user.email = "test@example.com"
    
    from amprenta_rag.api.dependencies import get_current_user
    app.dependency_overrides[get_current_user] = lambda: mock_user
    
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


class TestGetDeletionImpact:
    """Tests for GET /lifecycle/impact endpoint."""
    
    def test_rejects_invalid_entity_type(self, client):
        """Invalid entity type returns 400."""
        response = client.get(f"/api/v1/lifecycle/impact/invalid_type/{uuid4()}")
        assert response.status_code == 400
        assert "Invalid entity_type" in response.json()["detail"]
    
    @patch("amprenta_rag.api.routers.lifecycle.calculate_deletion_impact")
    def test_returns_impact_for_valid_entity(self, mock_impact, client):
        """Valid entity returns impact data."""
        entity_id = uuid4()
        mock_impact.return_value = {
            "entity": {"type": "dataset", "id": str(entity_id), "name": "Test"},
            "impact": {"features": 10},
            "blocking_references": [],
            "can_delete": True,
        }
        
        response = client.get(f"/api/v1/lifecycle/impact/dataset/{entity_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["can_delete"] is True


class TestUpdateStatus:
    """Tests for POST /lifecycle/status endpoint."""
    
    def test_rejects_invalid_entity_type(self, client):
        """Invalid entity type returns 400."""
        response = client.post("/api/v1/lifecycle/status", json={
            "entity_type": "invalid",
            "entity_id": str(uuid4()),
            "new_status": "active"
        })
        assert response.status_code == 400
    
    def test_rejects_invalid_status(self, client):
        """Invalid status returns 400."""
        response = client.post("/api/v1/lifecycle/status", json={
            "entity_type": "dataset",
            "entity_id": str(uuid4()),
            "new_status": "invalid_status"
        })
        assert response.status_code == 400
    
    @patch("amprenta_rag.api.routers.lifecycle.update_lifecycle_status")
    def test_successful_update(self, mock_update, client):
        """Successful update returns success response."""
        mock_update.return_value = (True, "Status updated")
        
        response = client.post("/api/v1/lifecycle/status", json={
            "entity_type": "dataset",
            "entity_id": str(uuid4()),
            "new_status": "quarantined",
            "reason": "Data quality review"
        })
        
        assert response.status_code == 200
        assert response.json()["success"] is True


class TestBulkOperations:
    """Tests for bulk operation endpoints."""
    
    def test_bulk_status_rejects_over_100_entities(self, client):
        """Bulk operations reject >100 entities."""
        entity_ids = [str(uuid4()) for _ in range(101)]
        
        response = client.post("/api/v1/lifecycle/bulk/status", json={
            "entity_type": "dataset",
            "entity_ids": entity_ids,
            "new_status": "archived"
        })
        
        assert response.status_code == 400
        assert "Maximum 100" in response.json()["detail"]
    
    def test_bulk_archive_requires_confirmation(self, client):
        """Bulk archive requires confirmed=true."""
        response = client.post("/api/v1/lifecycle/bulk/archive", json={
            "entity_type": "dataset",
            "entity_ids": [str(uuid4())],
            "reason": "Cleanup",
            "confirmed": False
        })
        
        assert response.status_code == 400
        assert "confirmed=true" in response.json()["detail"]
    
    @patch("amprenta_rag.api.routers.lifecycle.bulk_delete_preview")
    def test_bulk_preview_returns_aggregated_impact(self, mock_preview, client):
        """Bulk preview returns aggregated impact."""
        mock_preview.return_value = {
            "entity_type": "dataset",
            "entity_count": 2,
            "total_impact": {"features": 30},
            "blocking_entities": [],
            "can_proceed": True,
        }
        
        response = client.post("/api/v1/lifecycle/bulk/preview", json={
            "entity_type": "dataset",
            "entity_ids": [str(uuid4()), str(uuid4())]
        })
        
        assert response.status_code == 200
        assert response.json()["can_proceed"] is True
