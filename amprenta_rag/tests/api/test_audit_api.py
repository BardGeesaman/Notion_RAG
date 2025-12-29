"""Tests for audit API endpoints."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestGetAuditTrail:
    """Tests for GET /api/v1/audit/{entity_type}/{entity_id} endpoint."""

    @patch("amprenta_rag.api.routers.audit.get_audit_trail")
    def test_get_audit_trail_success(self, mock_get_trail):
        """Test successful audit trail retrieval."""
        mock_log = MagicMock()
        mock_log.entity_type = "dataset"
        mock_log.entity_id = "test-id"
        mock_log.action = "update"
        mock_log.username = "testuser"
        mock_log.timestamp = datetime(2024, 1, 1)
        mock_log.old_checksum = "abc123"
        mock_log.new_checksum = "def456"
        mock_log.changes = {"name": {"old": "old", "new": "new"}}
        
        mock_get_trail.return_value = [mock_log]
        
        response = client.get("/api/v1/audit/dataset/test-id")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["action"] == "update"


class TestVerifyIntegrity:
    """Tests for POST /api/v1/audit/{entity_type}/{entity_id}/verify endpoint."""

    @patch("amprenta_rag.api.routers.audit.verify_integrity")
    def test_verify_integrity_success(self, mock_verify):
        """Test successful integrity verification."""
        mock_verify.return_value = True
        
        response = client.post(
            "/api/v1/audit/dataset/test-id/verify",
            json={"current_data": {"name": "test"}},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["verified"] is True

    @patch("amprenta_rag.api.routers.audit.verify_integrity")
    def test_verify_integrity_failed(self, mock_verify):
        """Test failed integrity verification."""
        mock_verify.return_value = False
        
        response = client.post(
            "/api/v1/audit/dataset/test-id/verify",
            json={"current_data": {"name": "modified"}},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["verified"] is False


class TestGetRecentChanges:
    """Tests for GET /api/v1/audit/recent endpoint."""

    def test_get_recent_changes(self):
        """Test recent changes retrieval."""
        mock_log = MagicMock()
        mock_log.entity_type = "experiment"
        mock_log.entity_id = "exp-id"
        mock_log.action = "create"
        mock_log.username = "scientist"
        mock_log.timestamp = datetime(2024, 1, 1)
        mock_log.old_checksum = None
        mock_log.new_checksum = "xyz789"
        mock_log.changes = None
        
        mock_db = MagicMock()
        mock_db.query.return_value.order_by.return_value.offset.return_value.limit.return_value.all.return_value = [mock_log]
        
        def mock_get_db():
            yield mock_db
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get("/api/v1/audit/recent?limit=10")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
        finally:
            app.dependency_overrides.clear()


class TestGetUserChanges:
    """Tests for GET /api/v1/audit/user/{user_id} endpoint."""

    @patch("amprenta_rag.api.routers.audit.AuditLog")
    def test_get_user_changes(self, mock_audit_log_class):
        """Test user changes retrieval."""
        mock_log = MagicMock()
        mock_log.entity_type = "compound"
        mock_log.entity_id = "cmpd-id"
        mock_log.action = "delete"
        mock_log.username = "admin"
        mock_log.timestamp = datetime(2024, 1, 1)
        mock_log.old_checksum = "old123"
        mock_log.new_checksum = None
        mock_log.changes = None
        
        mock_db = MagicMock()
        
        # Setup query chain
        mock_query = MagicMock()
        mock_query.filter.return_value = mock_query
        mock_query.order_by.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = [mock_log]
        
        mock_db.query.return_value = mock_query
        
        def mock_get_db():
            yield mock_db
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get("/api/v1/audit/user/user-123?limit=50")
            
            assert response.status_code == 200
            data = response.json()
            assert isinstance(data, list)
        finally:
            app.dependency_overrides.clear()

