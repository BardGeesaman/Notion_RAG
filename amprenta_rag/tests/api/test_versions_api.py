"""Tests for version management API endpoints."""

import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch
from datetime import datetime

from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_database_session


class TestVersionsAPI:
    """Test cases for version management API endpoints."""

    @pytest.fixture
    def mock_user(self):
        """Mock authenticated user."""
        user = MagicMock()
        user.id = uuid4()
        user.role = "user"
        return user

    @pytest.fixture
    def client(self, mock_user):
        """Test client with mocked authentication."""
        def override_get_current_user():
            return mock_user

        def override_get_database_session():
            return MagicMock()

        app.dependency_overrides[get_current_user] = override_get_current_user
        app.dependency_overrides[get_database_session] = override_get_database_session
        
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    def test_list_versions_success(self, client):
        """Test successful listing of versions for an entity."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        # Mock version objects
        mock_version1 = MagicMock()
        mock_version1.id = uuid4()
        mock_version1.entity_type = entity_type
        mock_version1.entity_id = entity_id
        mock_version1.version_number = 2
        mock_version1.checksum_sha256 = "abc123"
        mock_version1.created_by = uuid4()
        mock_version1.created_at = datetime.now()
        mock_version1.change_summary = "Updated dataset"
        
        mock_version2 = MagicMock()
        mock_version2.id = uuid4()
        mock_version2.entity_type = entity_type
        mock_version2.entity_id = entity_id
        mock_version2.version_number = 1
        mock_version2.checksum_sha256 = "def456"
        mock_version2.created_by = uuid4()
        mock_version2.created_at = datetime.now()
        mock_version2.change_summary = "Initial version"
        
        with patch("amprenta_rag.api.routers.versions.get_versions") as mock_get_versions:
            mock_get_versions.return_value = [mock_version1, mock_version2]
            
            response = client.get(f"/api/v1/versions/{entity_type}/{entity_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 2
            assert data[0]["version_number"] == 2
            assert data[1]["version_number"] == 1
            mock_get_versions.assert_called_once()

    def test_list_versions_empty(self, client):
        """Test listing versions when no versions exist."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        with patch("amprenta_rag.api.routers.versions.get_versions") as mock_get_versions:
            mock_get_versions.return_value = []
            
            response = client.get(f"/api/v1/versions/{entity_type}/{entity_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert data == []

    def test_list_versions_invalid_entity_type(self, client):
        """Test listing versions with invalid entity type."""
        entity_id = uuid4()
        entity_type = "invalid_type"
        
        response = client.get(f"/api/v1/versions/{entity_type}/{entity_id}")
        
        assert response.status_code == 400
        assert "not versionable" in response.json()["detail"]

    def test_get_version_detail_success(self, client):
        """Test successful retrieval of version detail."""
        version_id = uuid4()
        
        mock_version = MagicMock()
        mock_version.id = version_id
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 1
        mock_version.checksum_sha256 = "abc123"
        mock_version.created_by = uuid4()
        mock_version.created_at = datetime.now()
        mock_version.change_summary = "Initial version"
        mock_version.data_snapshot = {"id": "123", "name": "Test Dataset"}
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version:
            mock_get_version.return_value = mock_version
            
            response = client.get(f"/api/v1/versions/{version_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == str(version_id)
            assert data["version_number"] == 1
            assert data["data_snapshot"]["name"] == "Test Dataset"
            mock_get_version.assert_called_once()

    def test_get_version_detail_not_found(self, client):
        """Test version detail retrieval when version doesn't exist."""
        version_id = uuid4()
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version:
            mock_get_version.return_value = None
            
            response = client.get(f"/api/v1/versions/{version_id}")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]

    def test_create_version_success(self, client, mock_user):
        """Test successful creation of a manual version."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        mock_version = MagicMock()
        mock_version.id = uuid4()
        mock_version.entity_type = entity_type
        mock_version.entity_id = entity_id
        mock_version.version_number = 1
        mock_version.checksum_sha256 = "abc123"
        mock_version.created_by = mock_user.id
        mock_version.created_at = datetime.now()
        mock_version.change_summary = "Manual snapshot"
        
        request_data = {
            "data": {"id": str(entity_id), "name": "Test Dataset"},
            "change_summary": "Manual snapshot"
        }
        
        with patch("amprenta_rag.api.routers.versions.create_version") as mock_create_version:
            mock_create_version.return_value = mock_version
            
            response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}", json=request_data)
            
            assert response.status_code == 201
            data = response.json()
            assert data["version_number"] == 1
            assert data["change_summary"] == "Manual snapshot"
            mock_create_version.assert_called_once()

    def test_create_version_invalid_entity_type(self, client):
        """Test version creation with invalid entity type."""
        entity_id = uuid4()
        entity_type = "invalid_type"
        
        request_data = {
            "data": {"id": str(entity_id), "name": "Test Dataset"},
            "change_summary": "Manual snapshot"
        }
        
        response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}", json=request_data)
        
        assert response.status_code == 400
        assert "not versionable" in response.json()["detail"]

    def test_create_version_service_error(self, client):
        """Test version creation when service raises an error."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        request_data = {
            "data": {"id": str(entity_id), "name": "Test Dataset"},
            "change_summary": "Manual snapshot"
        }
        
        with patch("amprenta_rag.api.routers.versions.create_version") as mock_create_version:
            mock_create_version.side_effect = ValueError("Invalid data")
            
            response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}", json=request_data)
            
            assert response.status_code == 400
            assert "Invalid data" in response.json()["detail"]

    def test_compare_versions_success(self, client):
        """Test successful version comparison."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        # Mock version objects
        mock_version1 = MagicMock()
        mock_version1.version_number = 1
        mock_version1.data_snapshot = {"name": "Original Dataset", "size": 100}
        
        mock_version2 = MagicMock()
        mock_version2.version_number = 2
        mock_version2.data_snapshot = {"name": "Updated Dataset", "size": 150, "description": "New desc"}
        
        request_data = {
            "version_a": 1,
            "version_b": 2
        }
        
        with patch("amprenta_rag.api.routers.versions.get_versions") as mock_get_versions, \
             patch("amprenta_rag.api.routers.versions.compare_versions") as mock_compare_versions:
            
            mock_get_versions.return_value = [mock_version2, mock_version1]  # Newest first
            mock_compare_versions.return_value = {
                "added": {"description": "New desc"},
                "removed": {},
                "changed": {"name": {"old": "Original Dataset", "new": "Updated Dataset"}, "size": {"old": 100, "new": 150}}
            }
            
            response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}/compare", json=request_data)
            
            assert response.status_code == 200
            data = response.json()
            assert data["version_a"] == 1
            assert data["version_b"] == 2
            assert "description" in data["added"]
            assert "name" in data["changed"]
            mock_compare_versions.assert_called_once()

    def test_compare_versions_missing_version(self, client):
        """Test version comparison when one version doesn't exist."""
        entity_id = uuid4()
        entity_type = "dataset"
        
        mock_version1 = MagicMock()
        mock_version1.version_number = 1
        mock_version1.data_snapshot = {"name": "Original Dataset"}
        
        request_data = {
            "version_a": 1,
            "version_b": 999  # Non-existent version
        }
        
        with patch("amprenta_rag.api.routers.versions.get_versions") as mock_get_versions:
            mock_get_versions.return_value = [mock_version1]
            
            response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}/compare", json=request_data)
            
            assert response.status_code == 404
            assert "Version 999 not found" in response.json()["detail"]

    def test_compare_versions_invalid_entity_type(self, client):
        """Test version comparison with invalid entity type."""
        entity_id = uuid4()
        entity_type = "invalid_type"
        
        request_data = {
            "version_a": 1,
            "version_b": 2
        }
        
        response = client.post(f"/api/v1/versions/{entity_type}/{entity_id}/compare", json=request_data)
        
        assert response.status_code == 400
        assert "not versionable" in response.json()["detail"]

    def test_restore_version_success(self, client, mock_user):
        """Test successful version restore by admin with confirmation."""
        version_id = uuid4()
        mock_user.role = "admin"  # Make user admin
        mock_user.username = "admin_user"
        
        # Mock the version to restore
        mock_version = MagicMock()
        mock_version.id = version_id
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 2
        mock_version.data_snapshot = {"name": "Previous Dataset", "size": 100}
        
        # Mock the new version created by rollback
        mock_new_version = MagicMock()
        mock_new_version.version_number = 4
        mock_new_version.change_summary = "Rolled back to version 2: Rolling back problematic changes"
        
        request_data = {
            "confirm": True,
            "reason": "Rolling back problematic changes"
        }
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version, \
             patch("amprenta_rag.api.routers.versions.rollback_to_version") as mock_rollback, \
             patch("amprenta_rag.auth.audit.log_action") as mock_log_action:
            
            mock_get_version.return_value = mock_version
            mock_rollback.return_value = mock_new_version
            
            response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
            
            assert response.status_code == 200
            data = response.json()
            assert data["success"] is True
            assert data["restored_from_version"] == 2
            assert data["new_version_number"] == 4
            assert "Rolling back problematic changes" in data["message"]
            
            # Verify audit logging was called
            mock_log_action.assert_called_once()
            audit_call = mock_log_action.call_args[1]
            assert audit_call["action"] == "restore_version"
            assert audit_call["details"]["restored_from_version"] == 2
            assert audit_call["details"]["new_version_number"] == 4
            assert audit_call["details"]["reason"] == "Rolling back problematic changes"

    def test_restore_version_no_confirm(self, client, mock_user):
        """Test restore fails without confirm=True."""
        version_id = uuid4()
        mock_user.role = "admin"
        
        request_data = {
            "confirm": False,  # Missing confirmation
            "reason": "Test restore"
        }
        
        response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
        
        assert response.status_code == 400
        assert "confirm=True" in response.json()["detail"]

    def test_restore_version_non_admin(self, client, mock_user):
        """Test restore fails for non-admin user."""
        version_id = uuid4()
        mock_user.role = "researcher"  # Non-admin role
        
        request_data = {
            "confirm": True,
            "reason": "Test restore"
        }
        
        response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
        
        assert response.status_code == 403
        assert "Admin access required" in response.json()["detail"]

    def test_restore_version_not_found(self, client, mock_user):
        """Test restore fails when version doesn't exist."""
        version_id = uuid4()
        mock_user.role = "admin"
        
        request_data = {
            "confirm": True,
            "reason": "Test restore"
        }
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version:
            mock_get_version.return_value = None
            
            response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]

    def test_restore_version_creates_audit(self, client, mock_user):
        """Test that restore operation creates proper audit log entry."""
        version_id = uuid4()
        mock_user.role = "admin"
        mock_user.username = "test_admin"
        
        mock_version = MagicMock()
        mock_version.id = version_id
        mock_version.entity_type = "experiment"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 1
        
        mock_new_version = MagicMock()
        mock_new_version.version_number = 3
        mock_new_version.change_summary = "Rolled back to version 1: Emergency rollback"
        
        request_data = {
            "confirm": True,
            "reason": "Emergency rollback"
        }
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version, \
             patch("amprenta_rag.api.routers.versions.rollback_to_version") as mock_rollback, \
             patch("amprenta_rag.auth.audit.log_action") as mock_log_action:
            
            mock_get_version.return_value = mock_version
            mock_rollback.return_value = mock_new_version
            
            response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
            
            assert response.status_code == 200
            
            # Verify audit log was called with correct parameters
            mock_log_action.assert_called_once()
            audit_call = mock_log_action.call_args
            
            # Check positional and keyword arguments
            assert "restore_version" in str(audit_call)
            assert str(mock_user.id) in str(audit_call)
            assert "test_admin" in str(audit_call)
            assert "experiment" in str(audit_call)
            assert str(mock_version.entity_id) in str(audit_call)

    def test_restore_version_with_reason(self, client, mock_user):
        """Test restore includes reason in change summary and response."""
        version_id = uuid4()
        mock_user.role = "admin"
        mock_user.username = "admin_user"
        
        mock_version = MagicMock()
        mock_version.id = version_id
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 1
        
        mock_new_version = MagicMock()
        mock_new_version.version_number = 2
        mock_new_version.change_summary = "Rolled back to version 1: Data corruption detected in current version"
        
        request_data = {
            "confirm": True,
            "reason": "Data corruption detected in current version"
        }
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version, \
             patch("amprenta_rag.api.routers.versions.rollback_to_version") as mock_rollback, \
             patch("amprenta_rag.auth.audit.log_action") as mock_log_action:
            
            mock_get_version.return_value = mock_version
            mock_rollback.return_value = mock_new_version
            
            response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
            
            assert response.status_code == 200
            data = response.json()
            
            # Verify reason is included in message
            assert "Rolled back to version 1" in data["message"]
            assert "Data corruption detected in current version" in data["message"]
            
            # Verify rollback was called with reason
            mock_rollback.assert_called_once()
            rollback_call = mock_rollback.call_args[1]
            assert rollback_call["reason"] == "Data corruption detected in current version"

    def test_restore_version_service_error(self, client, mock_user):
        """Test restore handles service errors gracefully."""
        version_id = uuid4()
        mock_user.role = "admin"
        
        mock_version = MagicMock()
        mock_version.id = version_id
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 1
        
        request_data = {
            "confirm": True,
            "reason": "Test restore"
        }
        
        with patch("amprenta_rag.api.routers.versions.get_version") as mock_get_version, \
             patch("amprenta_rag.api.routers.versions.rollback_to_version") as mock_rollback:
            
            mock_get_version.return_value = mock_version
            mock_rollback.side_effect = Exception("Database error")
            
            response = client.post(f"/api/v1/versions/restore/{version_id}", json=request_data)
            
            assert response.status_code == 500
            assert "Failed to restore version" in response.json()["detail"]
