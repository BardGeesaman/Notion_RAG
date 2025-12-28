"""Tests for sharing service functionality."""

import pytest
from unittest.mock import MagicMock
from uuid import uuid4

from amprenta_rag.models.auth import EntityShare
from amprenta_rag.services.sharing import (
    share_entity,
    unshare_entity,
    check_share_permission,
    list_shares_for_entity,
)


@pytest.fixture
def mock_db_session():
    """Mock database session."""
    session = MagicMock()
    # Make the session act as its own context manager
    session.__enter__ = MagicMock(return_value=session)
    session.__exit__ = MagicMock(return_value=None)
    return session


def test_share_entity_with_user(mock_db_session):
    """Test sharing entity with a user."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    user_id = uuid4()
    shared_by_id = uuid4()
    permission = "edit"
    
    # Configure mock query chain for checking existing share
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = None  # No existing share
    mock_db_session.query.return_value = mock_query
    
    # Act
    result = share_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        user_id=user_id,
        permission=permission,
        shared_by=shared_by_id,
        db=mock_db_session
    )
    
    # Assert
    assert result is not None
    mock_db_session.add.assert_called_once()
    mock_db_session.commit.assert_called_once()
    mock_db_session.expunge.assert_called_once()


def test_share_entity_with_team(mock_db_session):
    """Test sharing entity with a team."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    team_id = uuid4()
    shared_by_id = uuid4()
    permission = "view"
    
    # Configure mock query chain for checking existing share
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = None  # No existing share
    mock_db_session.query.return_value = mock_query
    
    # Act
    result = share_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        team_id=team_id,
        permission=permission,
        shared_by=shared_by_id,
        db=mock_db_session
    )
    
    # Assert
    assert result is not None
    mock_db_session.add.assert_called_once()
    mock_db_session.commit.assert_called_once()


def test_unshare_entity(mock_db_session):
    """Test unsharing an entity."""
    # Arrange
    entity_type = "compound"
    entity_id = uuid4()
    user_id = uuid4()
    shared_by_id = uuid4()
    
    # Mock existing share
    mock_share = EntityShare(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        shared_with_user_id=user_id,
        shared_with_team_id=None,
        permission="view",
        shared_by_id=shared_by_id
    )
    
    # Configure mock query chain for finding existing share
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = mock_share  # Return the share to delete
    mock_db_session.query.return_value = mock_query
    
    # Act
    result = unshare_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        user_id=user_id,
        db=mock_db_session
    )
    
    # Assert
    assert result is True
    mock_db_session.delete.assert_called_once_with(mock_share)
    mock_db_session.commit.assert_called_once()


def test_check_share_permission(mock_db_session):
    """Test checking share permissions."""
    # Arrange
    user_id = uuid4()
    entity_type = "signature"
    entity_id = uuid4()
    required_permission = "edit"
    
    # Mock existing share with edit permission
    mock_share = EntityShare(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        shared_with_user_id=user_id,
        shared_with_team_id=None,
        permission="edit",
        shared_by_id=uuid4()
    )
    
    # Configure mock query chains
    # First query for team membership (TeamMember)
    mock_team_query = MagicMock()
    mock_team_filter = MagicMock()
    mock_team_query.filter.return_value = mock_team_filter
    mock_team_filter.subquery.return_value = MagicMock()
    
    # Second query for entity shares (EntityShare)
    mock_share_query = MagicMock()
    mock_share_filter = MagicMock()
    mock_share_query.filter.return_value = mock_share_filter
    mock_share_filter.all.return_value = [mock_share]
    
    # Configure query to return different mocks based on call order
    mock_db_session.query.side_effect = [mock_team_query, mock_share_query]
    
    # Act
    result = check_share_permission(
        user_id=user_id,
        entity_type=entity_type,
        entity_id=entity_id,
        required_permission=required_permission,
        db=mock_db_session
    )
    
    # Assert
    assert result is True


def test_check_share_permission_insufficient(mock_db_session):
    """Test checking share permissions with insufficient access."""
    # Arrange
    user_id = uuid4()
    entity_type = "dataset"
    entity_id = uuid4()
    required_permission = "admin"
    
    # Mock team membership
    team_ids_subquery = MagicMock()
    mock_db_session.query.return_value.filter.return_value.subquery.return_value = team_ids_subquery
    
    # Mock existing share with only view permission
    mock_share = EntityShare(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        shared_with_user_id=user_id,
        shared_with_team_id=None,
        permission="view",
        shared_by_id=uuid4()
    )
    
    mock_db_session.query.return_value.filter.return_value.all.return_value = [mock_share]
    
    # Act
    result = check_share_permission(
        user_id=user_id,
        entity_type=entity_type,
        entity_id=entity_id,
        required_permission=required_permission,
        db=mock_db_session
    )
    
    # Assert
    assert result is False


def test_list_shares_for_entity(mock_db_session):
    """Test listing shares for an entity."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    
    mock_shares = [
        EntityShare(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            shared_with_user_id=uuid4(),
            shared_with_team_id=None,
            permission="view",
            shared_by_id=uuid4()
        ),
        EntityShare(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            shared_with_user_id=None,
            shared_with_team_id=uuid4(),
            permission="edit",
            shared_by_id=uuid4()
        )
    ]
    
    # Configure mock query chain for listing shares
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.all.return_value = mock_shares
    mock_db_session.query.return_value = mock_query
    
    # Act
    result = list_shares_for_entity(
        entity_type=entity_type,
        entity_id=entity_id,
        db=mock_db_session
    )
    
    # Assert
    assert len(result) == 2
    assert all(share.entity_type == entity_type for share in mock_shares)
    assert all(share.entity_id == entity_id for share in mock_shares)
    
    # Verify expunge was called for each share
    assert mock_db_session.expunge.call_count == 2
