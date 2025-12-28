"""Tests for entity reviews service functionality."""

import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.models.auth import EntityReview, User
from amprenta_rag.services.entity_reviews import (
    create_review,
    submit_for_review,
    assign_reviewer,
    decide_review,
)


@pytest.fixture
def mock_db_session():
    """Mock database session."""
    session = MagicMock()
    # Make the session act as its own context manager
    session.__enter__ = MagicMock(return_value=session)
    session.__exit__ = MagicMock(return_value=None)
    return session


def test_create_review(mock_db_session):
    """Test creating a new review."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    submitted_by_id = uuid4()
    
    # Configure mock query chain for checking existing review
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = None  # No existing review
    mock_db_session.query.return_value = mock_query
    
    # Act
    result = create_review(
        entity_type=entity_type,
        entity_id=entity_id,
        submitted_by_id=submitted_by_id,
        db=mock_db_session
    )
    
    # Assert
    assert result is not None
    mock_db_session.add.assert_called_once()
    mock_db_session.commit.assert_called_once()
    mock_db_session.expunge.assert_called_once()


def test_create_review_duplicate_active(mock_db_session):
    """Test creating review when active review already exists."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    submitted_by_id = uuid4()
    
    # Mock existing active review
    existing_review = EntityReview(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        reviewer_id=uuid4(),
        status="in_review",
        comments="Existing review"
    )
    
    mock_db_session.query.return_value.filter.return_value.first.return_value = existing_review
    
    # Act
    result = create_review(
        entity_type=entity_type,
        entity_id=entity_id,
        submitted_by_id=submitted_by_id,
        db=mock_db_session
    )
    
    # Assert
    assert result is None
    mock_db_session.add.assert_not_called()


def test_submit_for_review(mock_db_session):
    """Test submitting a review for approval."""
    # Arrange
    review_id = uuid4()
    
    # Mock existing review
    mock_review = EntityReview(
        id=review_id,
        entity_type="compound",
        entity_id=uuid4(),
        reviewer_id=uuid4(),
        status="draft",
        comments="Review created"
    )
    
    # Configure mock query chain for finding review
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = mock_review
    mock_db_session.query.return_value = mock_query
    
    with patch('amprenta_rag.services.entity_reviews._is_valid_transition') as mock_valid:
        mock_valid.return_value = True
        
        # Act
        result = submit_for_review(review_id=review_id, db=mock_db_session)
        
        # Assert
        assert result is not None
        mock_db_session.commit.assert_called()
        mock_db_session.expunge.assert_called()


def test_assign_reviewer(mock_db_session):
    """Test assigning a reviewer to a review."""
    # Arrange
    review_id = uuid4()
    reviewer_id = uuid4()
    
    # Mock existing review
    mock_review = EntityReview(
        id=review_id,
        entity_type="signature",
        entity_id=uuid4(),
        reviewer_id=None,  # Start with no reviewer
        status="submitted",
        comments="Review submitted"
    )
    
    # Mock reviewer user
    mock_reviewer = User(
        id=reviewer_id,
        username=f"reviewer_{uuid4().hex[:8]}",
        email=f"reviewer_{uuid4().hex[:8]}@test.com",
        password_hash="hash"
    )
    
    # Configure mock queries for both review and reviewer lookups
    def mock_query_side_effect(model):
        if model.__name__ == "EntityReview":
            query = MagicMock()
            filter_obj = MagicMock()
            query.filter.return_value = filter_obj
            filter_obj.first.return_value = mock_review
            return query
        elif model.__name__ == "User":
            query = MagicMock()
            filter_obj = MagicMock()
            query.filter.return_value = filter_obj
            filter_obj.first.return_value = mock_reviewer
            return query
        return MagicMock()
    
    mock_db_session.query.side_effect = mock_query_side_effect
    
    with patch('amprenta_rag.services.entity_reviews._is_valid_transition') as mock_valid:
        mock_valid.return_value = True
        
        # Act
        result = assign_reviewer(
            review_id=review_id,
            reviewer_id=reviewer_id,
            db=mock_db_session
        )
        
        # Assert
        assert result is not None
        assert mock_review.reviewer_id == reviewer_id
        assert mock_review.status == "in_review"
        mock_db_session.commit.assert_called_once()
        mock_db_session.expunge.assert_called_once()


def test_decide_review_approved(mock_db_session):
    """Test approving a review."""
    # Arrange
    review_id = uuid4()
    decision = "approved"
    comment = "Looks good to me!"
    
    # Mock existing review
    mock_review = EntityReview(
        id=review_id,
        entity_type="dataset",
        entity_id=uuid4(),
        reviewer_id=uuid4(),
        status="in_review",
        comments="Review in progress"
    )
    
    # Configure mock query chain for finding review
    mock_query = MagicMock()
    mock_filter = MagicMock()
    mock_query.filter.return_value = mock_filter
    mock_filter.first.return_value = mock_review
    mock_db_session.query.return_value = mock_query
    
    with patch('amprenta_rag.services.entity_reviews._is_valid_transition') as mock_valid:
        mock_valid.return_value = True
        
        # Act
        result = decide_review(
            review_id=review_id,
            decision=decision,
            comment=comment,
            db=mock_db_session
        )
        
        # Assert
        assert result is not None
        mock_db_session.commit.assert_called()
        mock_db_session.expunge.assert_called()


def test_decide_review_invalid_decision(mock_db_session):
    """Test deciding review with invalid decision."""
    # Arrange
    review_id = uuid4()
    decision = "invalid_decision"
    comment = "This should fail"
    
    # Act
    result = decide_review(
        review_id=review_id,
        decision=decision,
        comment=comment,
        db=mock_db_session
    )
    
    # Assert
    assert result is None


def test_assign_reviewer_invalid_transition(mock_db_session):
    """Test assigning reviewer with invalid status transition."""
    # Arrange
    review_id = uuid4()
    reviewer_id = uuid4()
    
    # Mock existing review in wrong status
    mock_review = EntityReview(
        id=review_id,
        entity_type="experiment",
        entity_id=uuid4(),
        reviewer_id=uuid4(),
        status="approved",  # Can't assign reviewer to approved review
        comments="Already approved"
    )
    
    mock_db_session.query.return_value.filter.return_value.first.return_value = mock_review
    
    with patch('amprenta_rag.services.entity_reviews._is_valid_transition') as mock_valid:
        mock_valid.return_value = False
        
        # Act
        result = assign_reviewer(
            review_id=review_id,
            reviewer_id=reviewer_id,
            db=mock_db_session
        )
        
        # Assert
        assert result is None
        mock_db_session.commit.assert_not_called()


def test_assign_reviewer_not_found(mock_db_session):
    """Test assigning non-existent reviewer."""
    # Arrange
    review_id = uuid4()
    reviewer_id = uuid4()
    
    # Mock existing review
    mock_review = EntityReview(
        id=review_id,
        entity_type="compound",
        entity_id=uuid4(),
        reviewer_id=None,
        status="submitted",
        comments="Ready for review"
    )
    
    # Configure mock queries for review (found) and reviewer (not found)
    def mock_query_side_effect(model):
        if model.__name__ == "EntityReview":
            query = MagicMock()
            filter_obj = MagicMock()
            query.filter.return_value = filter_obj
            filter_obj.first.return_value = mock_review
            return query
        elif model.__name__ == "User":
            query = MagicMock()
            filter_obj = MagicMock()
            query.filter.return_value = filter_obj
            filter_obj.first.return_value = None  # Reviewer not found
            return query
        return MagicMock()
    
    mock_db_session.query.side_effect = mock_query_side_effect
    
    with patch('amprenta_rag.services.entity_reviews._is_valid_transition') as mock_valid:
        mock_valid.return_value = True
        
        # Act
        result = assign_reviewer(
            review_id=review_id,
            reviewer_id=reviewer_id,
            db=mock_db_session
        )
        
        # Assert
        assert result is None
        mock_db_session.commit.assert_not_called()
