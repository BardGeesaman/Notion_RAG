"""Tests for entity permissions functionality."""

import pytest
from unittest.mock import MagicMock, patch, ANY
from uuid import uuid4

from amprenta_rag.models.auth import User
from amprenta_rag.auth.permissions import (
    can_view_entity,
    can_edit_entity,
    can_share_entity,
    is_admin,
)


@pytest.fixture
def mock_db_session():
    """Mock database session."""
    session = MagicMock()
    # Make the session act as its own context manager
    session.__enter__ = MagicMock(return_value=session)
    session.__exit__ = MagicMock(return_value=None)
    return session


def test_can_view_entity_owner(mock_db_session):
    """Test that entity owner can view entity."""
    # Arrange
    user_id = uuid4()
    entity_type = "dataset"
    entity_id = uuid4()
    
    # Mock user (not admin)
    mock_user = User(
        id=user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Mock dataset owned by user
    mock_dataset = MagicMock()
    mock_dataset.created_by_id = user_id
    
    # Configure mock queries for User and Dataset lookups
    def mock_query_side_effect(model):
        query = MagicMock()
        filter_obj = MagicMock()
        query.filter.return_value = filter_obj
        
        if hasattr(model, '__name__') and model.__name__ == "User":
            filter_obj.first.return_value = mock_user
        else:
            # Dataset/entity query
            filter_obj.first.return_value = mock_dataset
        
        return query
    
    mock_db_session.query.side_effect = mock_query_side_effect
    
    with patch('amprenta_rag.auth.permissions.check_share_permission') as mock_check_share:
        mock_check_share.return_value = False  # No share permission needed
        
        # Act
        result = can_view_entity(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            db=mock_db_session
        )
        
        # Assert
        assert result is True


def test_can_edit_entity_with_share(mock_db_session):
    """Test that user with edit share can edit entity."""
    # Arrange
    user_id = uuid4()
    entity_type = "experiment"
    entity_id = uuid4()
    
    # Mock user (not admin)
    mock_user = User(
        id=user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Mock experiment not owned by user
    mock_experiment = MagicMock()
    mock_experiment.created_by_id = uuid4()  # Different user
    
    mock_db_session.query.return_value.filter.return_value.first.side_effect = [
        mock_user,  # User query
        mock_experiment  # Experiment query
    ]
    
    with patch('amprenta_rag.auth.permissions.check_share_permission') as mock_check_share:
        mock_check_share.return_value = True  # Has edit permission via share
        
        # Act
        result = can_edit_entity(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            db=mock_db_session
        )
        
        # Assert
        assert result is True
        mock_check_share.assert_called_with(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            required_permission="edit",
            db=ANY
        )


def test_permission_hierarchy(mock_db_session):
    """Test permission hierarchy (admin > owner > share)."""
    # Arrange
    admin_user_id = uuid4()
    entity_type = "compound"
    entity_id = uuid4()
    
    # Mock admin user
    mock_admin = User(
        id=admin_user_id,
        username=f"admin_{uuid4().hex[:8]}",
        email=f"admin_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="admin"
    )
    
    # Configure mock to always return admin user for User queries
    def mock_query_side_effect(model):
        query = MagicMock()
        filter_obj = MagicMock()
        query.filter.return_value = filter_obj
        filter_obj.first.return_value = mock_admin
        return query
    
    mock_db_session.query.side_effect = mock_query_side_effect
    
    # Act - Test all permission types for admin
    can_view = can_view_entity(
        user_id=admin_user_id,
        entity_type=entity_type,
        entity_id=entity_id,
        db=mock_db_session
    )
    
    can_edit = can_edit_entity(
        user_id=admin_user_id,
        entity_type=entity_type,
        entity_id=entity_id,
        db=mock_db_session
    )
    
    can_share = can_share_entity(
        user_id=admin_user_id,
        entity_type=entity_type,
        entity_id=entity_id,
        db=mock_db_session
    )
    
    # Assert - Admin should have all permissions
    assert can_view is True
    assert can_edit is True
    assert can_share is True


def test_orphan_handling(mock_db_session):
    """Test handling of entities without owners (orphaned entities)."""
    # Arrange
    user_id = uuid4()
    entity_type = "signature"
    entity_id = uuid4()
    
    # Mock user (not admin)
    mock_user = User(
        id=user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Mock signature with no owner (None created_by_id)
    mock_signature = MagicMock()
    mock_signature.created_by_id = None
    
    mock_db_session.query.return_value.filter.return_value.first.side_effect = [
        mock_user,  # User query
        mock_signature  # Signature query
    ]
    
    with patch('amprenta_rag.auth.permissions.check_share_permission') as mock_check_share:
        mock_check_share.return_value = False  # No share permission
        
        # Act
        result = can_edit_entity(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            db=mock_db_session
        )
        
        # Assert - Should not be able to edit orphaned entity without share
        assert result is False


def test_is_admin_function(mock_db_session):
    """Test is_admin helper function."""
    # Arrange
    admin_user_id = uuid4()
    regular_user_id = uuid4()
    
    # Mock admin user
    mock_admin = User(
        id=admin_user_id,
        username=f"admin_{uuid4().hex[:8]}",
        email=f"admin_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="admin"
    )
    
    # Mock regular user
    mock_user = User(
        id=regular_user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Configure mock query chain to return different users based on call order
    call_count = 0
    def mock_query_side_effect(model):
        nonlocal call_count
        query = MagicMock()
        filter_obj = MagicMock()
        query.filter.return_value = filter_obj
        
        if call_count == 0:
            filter_obj.first.return_value = mock_admin
        else:
            filter_obj.first.return_value = mock_user
        call_count += 1
        
        return query
    
    mock_db_session.query.side_effect = mock_query_side_effect
    
    # Act & Assert
    assert is_admin(admin_user_id, db=mock_db_session) is True
    assert is_admin(regular_user_id, db=mock_db_session) is False


def test_can_share_entity_requires_admin_permission(mock_db_session):
    """Test that sharing requires admin-level permission."""
    # Arrange
    user_id = uuid4()
    entity_type = "dataset"
    entity_id = uuid4()
    
    # Mock user (not admin)
    mock_user = User(
        id=user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Mock dataset not owned by user
    mock_dataset = MagicMock()
    mock_dataset.created_by_id = uuid4()  # Different user
    
    mock_db_session.query.return_value.filter.return_value.first.side_effect = [
        mock_user,  # User query
        mock_dataset  # Dataset query
    ]
    
    with patch('amprenta_rag.auth.permissions.check_share_permission') as mock_check_share:
        mock_check_share.return_value = False  # No admin share permission
        
        # Act
        result = can_share_entity(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            db=mock_db_session
        )
        
        # Assert - Should not be able to share without admin permission
        assert result is False
        mock_check_share.assert_called_with(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            required_permission="admin",
            db=ANY
        )


def test_compound_ownership_special_case(mock_db_session):
    """Test special case for compound ownership (admin-only edit)."""
    # Arrange
    user_id = uuid4()
    entity_type = "compound"
    entity_id = uuid4()
    
    # Mock user (not admin)
    mock_user = User(
        id=user_id,
        username=f"user_{uuid4().hex[:8]}",
        email=f"user_{uuid4().hex[:8]}@test.com",
        password_hash="hash",
        role="researcher"
    )
    
    # Mock compound (compounds don't have created_by_id in current implementation)
    mock_compound = MagicMock()
    
    mock_db_session.query.return_value.filter.return_value.first.side_effect = [
        mock_user,  # User query
        mock_compound  # Compound query
    ]
    
    with patch('amprenta_rag.auth.permissions.check_share_permission') as mock_check_share:
        mock_check_share.return_value = False  # No share permission
        
        # Act
        result = can_edit_entity(
            user_id=user_id,
            entity_type=entity_type,
            entity_id=entity_id,
            db=mock_db_session
        )
        
        # Assert - Regular user should not be able to edit compounds
        assert result is False
