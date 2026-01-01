"""Integration tests for collaboration API with real database."""

import pytest
from datetime import datetime
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.models.auth import EntityShare
from amprenta_rag.database.models import Program


@pytest.mark.integration
class TestCollaborationAPIIntegration:
    """Integration tests for collaboration API endpoints."""

    def test_list_sessions_success(self, integration_client, db_session, test_user, timed_request):
        """Test listing collaboration sessions from real database."""
        # Create real entity shares (collaboration sessions) in database
        shares = []
        for i in range(3):
            share = EntityShare(
                id=uuid4(),
                entity_type="program",
                entity_id=uuid4(),
                shared_by_id=test_user.id,
                shared_with_id=uuid4(),
                permission_level="read",
                created_at=datetime.utcnow()
            )
            shares.append(share)
        
        db_session.add_all(shares)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            "/api/v1/collaboration/sessions",
            "test_list_sessions_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert len(data["sessions"]) >= 3

    def test_list_sessions_empty(self, integration_client, db_session, timed_request):
        """Test listing sessions when none exist."""
        response, benchmark = timed_request(
            "GET",
            "/api/v1/collaboration/sessions",
            "test_list_sessions_empty"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "sessions" in data

    def test_get_session_success(self, integration_client, db_session, test_user, test_program, timed_request):
        """Test retrieving specific collaboration session from database."""
        # Create real entity share in database
        share = EntityShare(
            id=uuid4(),
            entity_type="program",
            entity_id=test_program.id,
            shared_by_id=test_user.id,
            shared_with_id=uuid4(),
            permission_level="write",
            created_at=datetime.utcnow()
        )
        db_session.add(share)
        db_session.commit()
        db_session.refresh(share)
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/collaboration/sessions/{share.id}",
            "test_get_session_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["id"] == str(share.id)
        assert data["entity_type"] == "program"
        assert data["permission_level"] == "write"

    def test_get_session_not_found(self, integration_client, timed_request):
        """Test retrieving non-existent session."""
        fake_session_id = uuid4()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/collaboration/sessions/{fake_session_id}",
            "test_get_session_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    def test_invite_to_session_success(self, integration_client, db_session, test_user, test_program, timed_request):
        """Test inviting user to collaboration session with real database."""
        # Create collaboration session
        session = EntityShare(
            id=uuid4(),
            entity_type="program",
            entity_id=test_program.id,
            shared_by_id=test_user.id,
            shared_with_id=uuid4(),
            permission_level="read",
            created_at=datetime.utcnow()
        )
        db_session.add(session)
        db_session.commit()
        
        # Create invitee user
        invitee = test_user.__class__(
            id=uuid4(),
            username=f"invitee_{uuid4().hex[:8]}",
            email=f"invitee_{uuid4().hex[:8]}@test.com",
            password_hash="hashedpassword",
            role="researcher"
        )
        db_session.add(invitee)
        db_session.commit()
        
        response, benchmark = timed_request(
            "POST",
            f"/api/v1/collaboration/sessions/{session.id}/invite",
            "test_invite_to_session_success",
            json={
                "user_email": invitee.email,
                "permission_level": "write"
            }
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["status"] == "invited"

    def test_invite_to_session_nonexistent_user(self, integration_client, db_session, test_user, test_program, timed_request):
        """Test inviting non-existent user returns error."""
        # Create collaboration session
        session = EntityShare(
            id=uuid4(),
            entity_type="program",
            entity_id=test_program.id,
            shared_by_id=test_user.id,
            shared_with_id=uuid4(),
            permission_level="read",
            created_at=datetime.utcnow()
        )
        db_session.add(session)
        db_session.commit()
        
        response, benchmark = timed_request(
            "POST",
            f"/api/v1/collaboration/sessions/{session.id}/invite",
            "test_invite_nonexistent_user",
            json={
                "user_email": "nonexistent@example.com",
                "permission_level": "read"
            }
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    @patch('amprenta_rag.api.routers.collaboration.get_yjs_server_url')
    def test_yjs_server_url_configuration(self, mock_yjs_url, integration_client, test_user, timed_request):
        """Test YJS server URL configuration with mocked external service."""
        # Mock external YJS server configuration
        mock_yjs_url.return_value = "ws://localhost:1234"
        
        response, benchmark = timed_request(
            "GET",
            "/api/v1/collaboration/config",
            "test_yjs_server_config"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "yjs_server_url" in data

    def test_entity_share_permissions(self, integration_client, db_session, test_user, test_program, timed_request):
        """Test entity sharing permission levels with real database."""
        # Create shares with different permission levels
        permissions = ["read", "write", "admin"]
        shares = []
        
        for perm in permissions:
            share = EntityShare(
                id=uuid4(),
                entity_type="program",
                entity_id=test_program.id,
                shared_by_id=test_user.id,
                shared_with_id=uuid4(),
                permission_level=perm,
                created_at=datetime.utcnow()
            )
            shares.append(share)
        
        db_session.add_all(shares)
        db_session.commit()
        
        # Verify permission levels stored correctly
        for share in shares:
            db_session.refresh(share)
            assert share.permission_level in permissions
            
        # Test retrieval includes permission levels
        response, benchmark = timed_request(
            "GET",
            "/api/v1/collaboration/sessions",
            "test_permission_levels"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
