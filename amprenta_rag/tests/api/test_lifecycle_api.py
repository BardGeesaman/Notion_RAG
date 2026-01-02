"""Tests for lifecycle management API endpoints."""

import pytest
from uuid import uuid4
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


# Mock user for auth
def mock_user():
    user = MagicMock()
    user.id = uuid4()
    user.email = "test@test.com"
    user.role = "admin"
    return user


@pytest.fixture
def client():
    app.dependency_overrides[get_current_user] = mock_user
    yield TestClient(app)
    app.dependency_overrides.clear()


class TestBulkDeleteAPI:
    """Tests for bulk delete endpoint."""
    
    def test_delete_requires_confirmation(self, client):
        """DELETE without confirmation should fail."""
        response = client.request(
            "DELETE",
            "/api/v1/lifecycle/bulk/delete",
            json={
                "entity_type": "dataset",
                "entity_ids": [str(uuid4())],
                "reason": "test",
                "confirmed": False,
                "dry_run": False
            }
        )
        assert response.status_code == 400
    
    def test_delete_dry_run_allowed(self, client):
        """DELETE with dry_run should work without confirmation."""
        response = client.request(
            "DELETE",
            "/api/v1/lifecycle/bulk/delete",
            json={
                "entity_type": "dataset",
                "entity_ids": [str(uuid4())],
                "reason": "test",
                "confirmed": False,
                "dry_run": True
            }
        )
        assert response.status_code == 200
    
    def test_delete_invalid_entity_type(self, client):
        """Invalid entity type should return 400."""
        response = client.request(
            "DELETE",
            "/api/v1/lifecycle/bulk/delete",
            json={
                "entity_type": "invalid",
                "entity_ids": [str(uuid4())],
                "reason": "test",
                "confirmed": True
            }
        )
        assert response.status_code == 400


class TestOrphansAPI:
    """Tests for orphans endpoint."""
    
    def test_get_orphan_stats(self, client):
        """GET /orphans should return stats."""
        response = client.get("/api/v1/lifecycle/orphans")
        assert response.status_code == 200
        data = response.json()
        assert "features" in data
        assert "signatures" in data


class TestExportAPI:
    """Tests for export endpoints."""
    
    def test_export_entity(self, client):
        """POST /export should return export data."""
        response = client.post(
            "/api/v1/lifecycle/export",
            json={
                "entity_type": "dataset",
                "entity_id": str(uuid4())
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert "entity_type" in data
    
    def test_export_download(self, client):
        """POST /export/download should return ZIP."""
        response = client.post(
            "/api/v1/lifecycle/export/download",
            json={
                "entity_type": "dataset",
                "entity_id": str(uuid4())
            }
        )
        assert response.status_code == 200
        assert response.headers.get("content-type") == "application/zip"


class TestStatusAPI:
    """Tests for status update endpoint."""
    
    def test_update_status_invalid(self, client):
        """Invalid status should return 400."""
        response = client.post(
            "/api/v1/lifecycle/status",
            json={
                "entity_type": "dataset",
                "entity_id": str(uuid4()),
                "new_status": "invalid_status"
            }
        )
        assert response.status_code == 400