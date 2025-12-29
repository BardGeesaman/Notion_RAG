"""Tests for share links API endpoints."""

from __future__ import annotations

from datetime import datetime, timedelta
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestCreateShareLink:
    """Tests for POST /api/v1/share-links endpoint."""

    @patch("amprenta_rag.api.routers.share_links.generate_share_link")
    def test_create_success(self, mock_generate):
        """Test successful share link creation."""
        mock_link = MagicMock()
        mock_link.id = uuid4()
        mock_link.token = "a" * 64
        mock_link.dashboard_path = "/voila/dashboard.ipynb"
        mock_link.expires_at = datetime.utcnow() + timedelta(hours=24)
        mock_link.max_views = 10
        mock_link.view_count = 0
        mock_link.is_active = True
        mock_link.permissions = "view"
        mock_link.created_at = datetime.utcnow()
        
        mock_generate.return_value = mock_link
        
        # Mock current_user dependency
        mock_user = MagicMock()
        mock_user.id = uuid4()
        mock_user.company_id = uuid4()
        
        def mock_get_user():
            return mock_user
        
        from amprenta_rag.api.dependencies import get_current_user
        
        app.dependency_overrides[get_current_user] = mock_get_user
        
        try:
            response = client.post(
                "/api/v1/share-links",
                json={
                    "dashboard_path": "/voila/dashboard.ipynb",
                    "context": {"experiment_id": "123"},
                    "expires_in_hours": 24,
                    "max_views": 10,
                    "permissions": "view",
                },
            )
            
            assert response.status_code == 201
            data = response.json()
            assert len(data["token"]) == 64
        finally:
            app.dependency_overrides.clear()


class TestListShareLinks:
    """Tests for GET /api/v1/share-links endpoint."""

    @patch("amprenta_rag.api.routers.share_links.get_user_share_links")
    def test_list_links(self, mock_get):
        """Test listing user's share links."""
        mock_link = MagicMock()
        mock_link.id = uuid4()
        mock_link.token = "test_token_123"
        mock_link.dashboard_path = "/voila/test.ipynb"
        mock_link.expires_at = datetime.utcnow()
        mock_link.max_views = None
        mock_link.view_count = 0
        mock_link.is_active = True
        mock_link.permissions = "view"
        mock_link.created_at = datetime.utcnow()
        
        mock_get.return_value = [mock_link]
        
        mock_user = MagicMock()
        mock_user.id = uuid4()
        
        def mock_get_user():
            return mock_user
        
        from amprenta_rag.api.dependencies import get_current_user
        
        app.dependency_overrides[get_current_user] = mock_get_user
        
        try:
            response = client.get("/api/v1/share-links")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
        finally:
            app.dependency_overrides.clear()


class TestRevokeShareLink:
    """Tests for DELETE /api/v1/share-links/{id} endpoint."""

    @patch("amprenta_rag.api.routers.share_links.revoke_share_link")
    def test_revoke(self, mock_revoke):
        """Test revoking share link."""
        mock_revoke.return_value = True
        
        mock_user = MagicMock()
        mock_user.id = uuid4()
        
        def mock_get_user():
            return mock_user
        
        from amprenta_rag.api.dependencies import get_current_user
        
        app.dependency_overrides[get_current_user] = mock_get_user
        
        try:
            link_id = uuid4()
            response = client.delete(f"/api/v1/share-links/{link_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert data["revoked"] is True
        finally:
            app.dependency_overrides.clear()


class TestValidateToken:
    """Tests for GET /api/v1/share-links/{token}/validate endpoint."""

    @patch("amprenta_rag.api.routers.share_links.validate_share_link")
    def test_validate_success(self, mock_validate):
        """Test successful token validation (PUBLIC endpoint)."""
        mock_validate.return_value = {
            "dashboard_path": "/voila/dashboard.ipynb",
            "context": {"experiment_id": "123"},
            "permissions": "view",
        }
        
        # NO auth mock - this is a public endpoint
        response = client.get("/api/v1/share-links/test_token_abc/validate")
        
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["voila_url"] is not None

    @patch("amprenta_rag.api.routers.share_links.validate_share_link")
    def test_validate_expired(self, mock_validate):
        """Test validation of expired token."""
        mock_validate.return_value = None  # Invalid/expired
        
        response = client.get("/api/v1/share-links/expired_token/validate")
        
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False


class TestGetStats:
    """Tests for GET /api/v1/share-links/{id}/stats endpoint."""

    @patch("amprenta_rag.api.routers.share_links.get_share_link_stats")
    def test_get_stats(self, mock_stats):
        """Test getting share link statistics."""
        mock_stats.return_value = {
            "view_count": 5,
            "max_views": 10,
            "last_accessed_at": datetime.utcnow().isoformat(),
            "time_remaining_hours": 12.5,
            "is_active": True,
        }
        
        mock_user = MagicMock()
        mock_user.id = uuid4()
        
        def mock_get_user():
            return mock_user
        
        from amprenta_rag.api.dependencies import get_current_user
        
        app.dependency_overrides[get_current_user] = mock_get_user
        
        try:
            link_id = uuid4()
            response = client.get(f"/api/v1/share-links/{link_id}/stats")
            
            assert response.status_code == 200
            data = response.json()
            assert data["view_count"] == 5
        finally:
            app.dependency_overrides.clear()

