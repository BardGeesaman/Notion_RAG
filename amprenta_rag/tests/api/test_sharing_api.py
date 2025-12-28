"""Tests for sharing API endpoints."""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock, patch, ANY
from uuid import uuid4
from datetime import datetime, timezone

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.models.auth import User, EntityShare


# Test constants
TEST_USER_ID = uuid4()


def mock_current_user():
    """Mock current user for authentication."""
    return User(
        id=TEST_USER_ID,
        username=f"test_user_{uuid4().hex[:8]}",
        email=f"test_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )


def mock_db():
    """Mock database session."""
    return MagicMock()


@pytest.fixture
def client():
    """Test client with mocked dependencies."""
    app.dependency_overrides[get_current_user] = mock_current_user
    app.dependency_overrides[get_db] = mock_db
    
    try:
        with TestClient(app) as test_client:
            yield test_client
    finally:
        app.dependency_overrides.clear()


def test_share_endpoint(client):
    """Test sharing an entity via API."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    user_id = uuid4()
    
    share_data = {
        "shared_with_user_id": str(user_id),
        "permission": "edit"
    }
    
    mock_share = EntityShare(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        shared_with_user_id=user_id,
        shared_with_team_id=None,
        permission="edit",
        shared_by_id=TEST_USER_ID,
        created_at=datetime.now(timezone.utc)
    )
    
    # Mock the authorization and share_entity functions
    with patch('amprenta_rag.api.routers.sharing.can_share_entity') as mock_can_share, \
         patch('amprenta_rag.api.routers.sharing.share_entity') as mock_share_func:
        
        mock_can_share.return_value = True  # User has permission to share
        mock_share_func.return_value = mock_share
        
        # Act
        response = client.post(
            f"/api/v1/sharing/entities/{entity_type}/{entity_id}/share",
            json=share_data
        )
        
        # Assert
        assert response.status_code == 201
        data = response.json()
        assert data["entity_type"] == entity_type
        assert data["entity_id"] == str(entity_id)
        assert data["shared_with_user_id"] == str(user_id)
        assert data["permission"] == "edit"


def test_unshare_endpoint(client):
    """Test unsharing an entity via API."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    user_id = uuid4()
    
    # Mock the authorization and unshare_entity functions
    with patch('amprenta_rag.api.routers.sharing.can_share_entity') as mock_can_share, \
         patch('amprenta_rag.api.routers.sharing.unshare_entity') as mock_unshare_func:
        
        mock_can_share.return_value = True  # User has permission to unshare
        mock_unshare_func.return_value = True
        
        # Act
        response = client.delete(
            f"/api/v1/sharing/entities/{entity_type}/{entity_id}/share",
            params={"user_id": str(user_id)}
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["message"] == "Share removed successfully"


def test_list_shares_endpoint(client):
    """Test listing shares for an entity via API."""
    # Arrange
    entity_type = "compound"
    entity_id = uuid4()
    
    mock_shares = [
        EntityShare(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            shared_with_user_id=uuid4(),
            shared_with_team_id=None,
            permission="view",
            shared_by_id=TEST_USER_ID,
            created_at=datetime.now(timezone.utc)
        ),
        EntityShare(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            shared_with_user_id=None,
            shared_with_team_id=uuid4(),
            permission="edit",
            shared_by_id=TEST_USER_ID,
            created_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the authorization and list_shares_for_entity functions
    with patch('amprenta_rag.api.routers.sharing.can_share_entity') as mock_can_share, \
         patch('amprenta_rag.api.routers.sharing.list_shares_for_entity') as mock_list_func:
        
        mock_can_share.return_value = True  # User has permission to view shares
        mock_list_func.return_value = mock_shares
        
        # Act
        response = client.get(
            f"/api/v1/sharing/entities/{entity_type}/{entity_id}/shares"
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["shares"]) == 2
        assert data["shares"][0]["entity_type"] == entity_type
        assert data["shares"][1]["entity_type"] == entity_type


def test_my_shares_endpoint(client):
    """Test getting user's shared entities via API."""
    # Arrange
    mock_shares = [
        EntityShare(
            id=uuid4(),
            entity_type="dataset",
            entity_id=uuid4(),
            shared_with_user_id=TEST_USER_ID,
            shared_with_team_id=None,
            permission="view",
            shared_by_id=uuid4(),
            created_at=datetime.now(timezone.utc)
        ),
        EntityShare(
            id=uuid4(),
            entity_type="experiment",
            entity_id=uuid4(),
            shared_with_user_id=None,
            shared_with_team_id=uuid4(),
            permission="edit",
            shared_by_id=uuid4(),
            created_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the get_my_shares function
    with patch('amprenta_rag.api.routers.sharing.get_my_shares') as mock_get_func:
        mock_get_func.return_value = mock_shares
        
        # Act
        response = client.get("/api/v1/sharing/my-shares")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["shares"]) == 2
        assert any(share["entity_type"] == "dataset" for share in data["shares"])
        assert any(share["entity_type"] == "experiment" for share in data["shares"])


def test_check_permission_endpoint(client):
    """Test checking entity permission via API."""
    # Arrange
    entity_type = "signature"
    entity_id = uuid4()
    permission = "edit"
    
    # Mock the check_share_permission function
    with patch('amprenta_rag.api.routers.sharing.check_share_permission') as mock_check_func:
        mock_check_func.return_value = True
        
        # Act
        response = client.get(
            f"/api/v1/sharing/entities/{entity_type}/{entity_id}/check-permission",
            params={"permission": permission}
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["entity_type"] == entity_type
        assert data["entity_id"] == str(entity_id)
        assert data["permission"] == permission
        assert data["has_permission"] is True
        assert data["user_id"] == str(TEST_USER_ID)


def test_share_endpoint_invalid_entity_type(client):
    """Test sharing with invalid entity type."""
    # Arrange
    entity_type = "invalid_type"
    entity_id = uuid4()
    
    share_data = {
        "shared_with_user_id": str(uuid4()),
        "permission": "view"
    }
    
    # Act
    response = client.post(
        f"/api/v1/sharing/entities/{entity_type}/{entity_id}/share",
        json=share_data
    )
    
    # Assert
    assert response.status_code == 400
    assert "Invalid entity_type" in response.json()["detail"]


def test_unshare_endpoint_missing_params(client):
    """Test unsharing without required parameters."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    
    # Act - No user_id or team_id provided
    response = client.delete(
        f"/api/v1/sharing/entities/{entity_type}/{entity_id}/share"
    )
    
    # Assert
    assert response.status_code == 400
    assert "Either user_id or team_id query parameter must be provided" in response.json()["detail"]


def test_unshare_endpoint_not_found(client):
    """Test unsharing when share doesn't exist."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    user_id = uuid4()
    
    # Mock the authorization and unshare_entity functions
    with patch('amprenta_rag.api.routers.sharing.can_share_entity') as mock_can_share, \
         patch('amprenta_rag.api.routers.sharing.unshare_entity') as mock_unshare_func:
        
        mock_can_share.return_value = True  # User has permission to unshare
        mock_unshare_func.return_value = False  # But share doesn't exist
        
        # Act
        response = client.delete(
            f"/api/v1/sharing/entities/{entity_type}/{entity_id}/share",
            params={"user_id": str(user_id)}
        )
        
        # Assert
        assert response.status_code == 404
        assert "Share not found or failed to remove" in response.json()["detail"]


def test_check_permission_endpoint_invalid_permission(client):
    """Test checking permission with invalid permission level."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    invalid_permission = "invalid_permission"
    
    # Act
    response = client.get(
        f"/api/v1/sharing/entities/{entity_type}/{entity_id}/check-permission",
        params={"permission": invalid_permission}
    )
    
    # Assert
    assert response.status_code == 400
    assert "Invalid permission" in response.json()["detail"]


def test_my_shares_endpoint_with_filter(client):
    """Test getting user's shared entities with entity type filter."""
    # Arrange
    entity_type_filter = "dataset"
    
    mock_shares = [
        EntityShare(
            id=uuid4(),
            entity_type="dataset",
            entity_id=uuid4(),
            shared_with_user_id=TEST_USER_ID,
            shared_with_team_id=None,
            permission="view",
            shared_by_id=uuid4(),
            created_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the get_my_shares function
    with patch('amprenta_rag.api.routers.sharing.get_my_shares') as mock_get_func:
        mock_get_func.return_value = mock_shares
        
        # Act
        response = client.get(
            "/api/v1/sharing/my-shares",
            params={"entity_type": entity_type_filter}
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["shares"]) == 1
        assert data["shares"][0]["entity_type"] == entity_type_filter
        
        # Verify the service was called with the filter
        mock_get_func.assert_called_once_with(
            user_id=TEST_USER_ID,
            entity_type=entity_type_filter,
            db=ANY
        )
