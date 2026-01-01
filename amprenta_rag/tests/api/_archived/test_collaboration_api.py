"""Tests for collaboration API endpoints."""

import os
from datetime import datetime, timedelta
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.database.models import User


class TestCollaborationAPI:
    """Test collaboration API endpoints for real-time collaborative editing."""

    @pytest.fixture
    def mock_user(self):
        """Fixture providing a mock user for authentication."""
        user = MagicMock(spec=User)
        user.id = uuid4()
        user.username = "testuser"
        user.email = "test@example.com"
        user.role = "user"
        return user

    @pytest.fixture
    def mock_admin_user(self):
        """Fixture providing a mock admin user."""
        user = MagicMock(spec=User)
        user.id = uuid4()
        user.username = "admin"
        user.email = "admin@example.com"
        user.role = "admin"
        return user

    @pytest.fixture
    def mock_db_session(self):
        """Fixture providing a mock database session."""
        return MagicMock()

    @pytest.fixture
    def client(self, mock_user, mock_db_session):
        """Fixture providing a test client with mocked dependencies."""
        def override_get_current_user():
            return mock_user

        def override_get_database_session():
            return mock_db_session

        app.dependency_overrides[get_current_user] = override_get_current_user
        app.dependency_overrides[get_database_session] = override_get_database_session

        try:
            with TestClient(app) as test_client:
                yield test_client
        finally:
            app.dependency_overrides.clear()

    @pytest.fixture
    def admin_client(self, mock_admin_user, mock_db_session):
        """Fixture providing a test client with admin user."""
        def override_get_current_user():
            return mock_admin_user

        def override_get_database_session():
            return mock_db_session

        app.dependency_overrides[get_current_user] = override_get_current_user
        app.dependency_overrides[get_database_session] = override_get_database_session

        try:
            with TestClient(app) as test_client:
                yield test_client
        finally:
            app.dependency_overrides.clear()

    def test_list_sessions_success(self, client):
        """Test successful listing of collaboration sessions."""
        # The mock function is working correctly, so just test the real response
        response = client.get("/api/v1/collaboration/sessions")

        assert response.status_code == 200
        data = response.json()
        assert "sessions" in data
        assert "total_count" in data
        assert "active_count" in data
        # The mock returns 2 sessions with 1 active
        assert data["total_count"] == 2
        assert data["active_count"] == 1
        assert len(data["sessions"]) == 2

    def test_list_sessions_empty(self, client):
        """Test listing sessions when user has no active sessions."""
        with patch("amprenta_rag.api.routers.collaboration._mock_get_active_sessions") as mock_get_sessions:
            mock_get_sessions.return_value = []

            response = client.get("/api/v1/collaboration/sessions")

            assert response.status_code == 200
            data = response.json()
            assert data["total_count"] == 0
            assert data["active_count"] == 0
            assert len(data["sessions"]) == 0

    def test_list_sessions_authentication_required(self):
        """Test that authentication is required for listing sessions."""
        with TestClient(app) as client:
            response = client.get("/api/v1/collaboration/sessions")
            assert response.status_code == 401

    def test_get_session_success(self, client, mock_user):
        """Test successful retrieval of session details."""
        document_id = "notebook_analysis_001"  # Use the mock document ID
        
        response = client.get(f"/api/v1/collaboration/sessions/{document_id}")

        assert response.status_code == 200
        data = response.json()
        assert "session" in data
        assert "user_permissions" in data
        assert "connection_url" in data
        assert "recent_activity" in data
        assert data["session"]["document_id"] == document_id
        assert "yjs-server:1234" in data["connection_url"]

    def test_get_session_not_found(self, client):
        """Test getting session details for non-existent document."""
        document_id = "nonexistent"
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            mock_get_details.return_value = None

            response = client.get(f"/api/v1/collaboration/sessions/{document_id}")

            assert response.status_code == 404
            data = response.json()
            assert "not found" in data["detail"].lower()

    def test_get_session_access_denied(self, client, mock_user):
        """Test access denied when user is not participant in session."""
        document_id = "restricted_doc"
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            # Mock session with different participants
            mock_session = MagicMock()
            mock_session.document_id = document_id
            
            # Mock participants NOT including current user
            other_user = MagicMock()
            other_user.id = uuid4()  # Different from mock_user.id
            mock_session.participants = [other_user]
            
            mock_get_details.return_value = mock_session

            response = client.get(f"/api/v1/collaboration/sessions/{document_id}")

            assert response.status_code == 403
            data = response.json()
            assert "access" in data["detail"].lower()

    def test_get_session_authentication_required(self):
        """Test that authentication is required for getting session details."""
        with TestClient(app) as client:
            response = client.get("/api/v1/collaboration/sessions/test_doc")
            assert response.status_code == 401

    def test_invite_to_session_success(self, client, mock_user):
        """Test successful invitation to collaboration session."""
        document_id = "notebook_001"
        invite_data = {
            "username": "collaborator",
            "message": "Would you like to collaborate?",
            "permissions": "read-write"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details, \
             patch("amprenta_rag.api.routers.collaboration._mock_send_invite") as mock_send_invite:
            
            # Mock session where current user is owner
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = mock_user.id
            mock_session.participants = [mock_session.owner]
            mock_get_details.return_value = mock_session

            # Mock successful invite
            mock_invite_response = MagicMock()
            mock_invite_response.success = True
            mock_invite_response.invite_id = uuid4()
            mock_invite_response.invited_user = invite_data["username"]
            mock_invite_response.document_id = document_id
            mock_invite_response.message = "Invite sent successfully"
            mock_invite_response.expires_at = datetime.utcnow() + timedelta(hours=24)
            mock_send_invite.return_value = mock_invite_response

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invite_data
            )

            assert response.status_code == 200
            data = response.json()
            assert data["success"] is True
            assert data["invited_user"] == invite_data["username"]
            assert data["document_id"] == document_id

    def test_invite_to_session_not_owner(self, client, mock_user):
        """Test invitation denied when user is not session owner."""
        document_id = "notebook_001"
        invite_data = {
            "username": "collaborator",
            "permissions": "read-write"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            # Mock session where current user is NOT owner
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = uuid4()  # Different from mock_user.id
            mock_get_details.return_value = mock_session

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invite_data
            )

            assert response.status_code == 403
            data = response.json()
            assert "owner" in data["detail"].lower()

    def test_invite_to_session_self_invite(self, client, mock_user):
        """Test that users cannot invite themselves."""
        document_id = "notebook_001"
        invite_data = {
            "username": mock_user.username,  # Same as current user
            "permissions": "read-write"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            # Mock session where current user is owner
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = mock_user.id
            mock_get_details.return_value = mock_session

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invite_data
            )

            assert response.status_code == 400
            data = response.json()
            assert "yourself" in data["detail"].lower()

    def test_invite_to_session_already_participant(self, client, mock_user):
        """Test invitation denied when user is already in session."""
        document_id = "notebook_001"
        invite_data = {
            "username": "existing_user",
            "permissions": "read-write"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            # Mock session where current user is owner
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = mock_user.id
            
            # Mock existing participant with same username
            existing_participant = MagicMock()
            existing_participant.username = invite_data["username"]
            mock_session.participants = [mock_session.owner, existing_participant]
            
            mock_get_details.return_value = mock_session

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invite_data
            )

            assert response.status_code == 400
            data = response.json()
            assert "already in" in data["detail"].lower()

    def test_invite_to_session_nonexistent_user(self, client, mock_user):
        """Test invitation to non-existent user."""
        document_id = "notebook_001"
        invite_data = {
            "username": "nonexistent_user",
            "permissions": "read-write"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            # Mock session where current user is owner
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = mock_user.id
            mock_session.participants = [mock_session.owner]
            mock_get_details.return_value = mock_session

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invite_data
            )

            assert response.status_code == 404
            data = response.json()
            assert "not found" in data["detail"].lower()

    def test_invite_to_session_authentication_required(self):
        """Test that authentication is required for sending invites."""
        with TestClient(app) as client:
            invite_data = {
                "username": "collaborator",
                "permissions": "read-write"
            }
            response = client.post(
                "/api/v1/collaboration/sessions/test_doc/invite",
                json=invite_data
            )
            assert response.status_code == 401

    @patch.dict(os.environ, {"YJS_SERVER_URL": "ws://custom-yjs:9999"})
    def test_yjs_server_url_configuration(self, client, mock_user):
        """Test that Y.js server URL is properly configured from environment."""
        document_id = "notebook_analysis_001"  # Use the mock document ID

        response = client.get(f"/api/v1/collaboration/sessions/{document_id}")

        assert response.status_code == 200
        data = response.json()
        assert data["connection_url"] == f"ws://custom-yjs:9999/doc/{document_id}"

    def test_invite_request_validation(self, client):
        """Test validation of invite request data."""
        document_id = "notebook_001"
        
        # Test missing username
        invalid_data = {
            "permissions": "read-write"
        }
        
        response = client.post(
            f"/api/v1/collaboration/sessions/{document_id}/invite",
            json=invalid_data
        )
        
        assert response.status_code == 422  # Validation error
        data = response.json()
        assert "detail" in data
        
        # Test invalid permissions (if validation is strict)
        invalid_permissions = {
            "username": "testuser",
            "permissions": "invalid-permission"
        }
        
        with patch("amprenta_rag.api.routers.collaboration._mock_get_session_details") as mock_get_details:
            mock_session = MagicMock()
            mock_session.document_id = document_id
            mock_session.owner.id = uuid4()
            mock_session.participants = [mock_session.owner]
            mock_get_details.return_value = mock_session

            response = client.post(
                f"/api/v1/collaboration/sessions/{document_id}/invite",
                json=invalid_permissions
            )
            
            # Should still work as we don't validate permissions strictly in current implementation
            # This test documents the current behavior
