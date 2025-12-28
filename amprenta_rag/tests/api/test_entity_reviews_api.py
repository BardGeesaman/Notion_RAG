"""Tests for entity reviews API endpoints."""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock, patch
from uuid import uuid4
from datetime import datetime, timezone

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.models.auth import User, EntityReview


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


def test_create_review_endpoint(client):
    """Test creating a review via API."""
    # Arrange
    entity_type = "dataset"
    entity_id = uuid4()
    
    mock_review = EntityReview(
        id=uuid4(),
        entity_type=entity_type,
        entity_id=entity_id,
        reviewer_id=TEST_USER_ID,
        status="draft",
        comments="Review created",
        reviewed_at=datetime.now(timezone.utc)
    )
    
    # Mock the create_review function
    with patch('amprenta_rag.api.routers.entity_reviews.create_review') as mock_create_func:
        mock_create_func.return_value = mock_review
        
        # Act
        response = client.post(
            f"/api/v1/reviews/entities/{entity_type}/{entity_id}/reviews"
        )
        
        # Assert
        assert response.status_code == 201
        data = response.json()
        assert data["entity_type"] == entity_type
        assert data["entity_id"] == str(entity_id)
        assert data["reviewer_id"] == str(TEST_USER_ID)
        assert data["status"] == "draft"


def test_list_reviews_endpoint(client):
    """Test listing reviews for an entity via API."""
    # Arrange
    entity_type = "experiment"
    entity_id = uuid4()
    
    mock_reviews = [
        EntityReview(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            reviewer_id=TEST_USER_ID,
            status="approved",
            comments="Approved",
            reviewed_at=datetime.now(timezone.utc)
        ),
        EntityReview(
            id=uuid4(),
            entity_type=entity_type,
            entity_id=entity_id,
            reviewer_id=uuid4(),
            status="draft",
            comments="In progress",
            reviewed_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the list_reviews_for_entity function
    with patch('amprenta_rag.api.routers.entity_reviews.list_reviews_for_entity') as mock_list_func:
        mock_list_func.return_value = mock_reviews
        
        # Act
        response = client.get(
            f"/api/v1/reviews/entities/{entity_type}/{entity_id}/reviews"
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["reviews"]) == 2
        assert data["reviews"][0]["entity_type"] == entity_type
        assert data["reviews"][1]["entity_type"] == entity_type


def test_submit_review_endpoint(client):
    """Test submitting a review via API."""
    # Arrange
    review_id = uuid4()
    
    mock_review = EntityReview(
        id=review_id,
        entity_type="compound",
        entity_id=uuid4(),
        reviewer_id=TEST_USER_ID,
        status="submitted",
        comments="Review submitted for approval",
        reviewed_at=datetime.now(timezone.utc)
    )
    
    # Mock the submit_for_review function
    with patch('amprenta_rag.api.routers.entity_reviews.submit_for_review') as mock_submit_func:
        mock_submit_func.return_value = mock_review
        
        # Act
        response = client.post(f"/api/v1/reviews/reviews/{review_id}/submit")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(review_id)
        assert data["status"] == "submitted"


def test_assign_reviewer_endpoint(client):
    """Test assigning a reviewer via API."""
    # Arrange
    review_id = uuid4()
    reviewer_id = uuid4()
    
    assignment_data = {
        "reviewer_id": str(reviewer_id)
    }
    
    mock_review = EntityReview(
        id=review_id,
        entity_type="signature",
        entity_id=uuid4(),
        reviewer_id=reviewer_id,
        status="in_review",
        comments="Review assigned to reviewer",
        reviewed_at=datetime.now(timezone.utc)
    )
    
    # Mock the assign_reviewer function
    with patch('amprenta_rag.api.routers.entity_reviews.assign_reviewer') as mock_assign_func:
        mock_assign_func.return_value = mock_review
        
        # Act
        response = client.post(
            f"/api/v1/reviews/reviews/{review_id}/assign",
            json=assignment_data
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(review_id)
        assert data["reviewer_id"] == str(reviewer_id)
        assert data["status"] == "in_review"


def test_decide_review_endpoint(client):
    """Test making a review decision via API."""
    # Arrange
    review_id = uuid4()
    
    decision_data = {
        "decision": "approved",
        "comment": "This looks great!"
    }
    
    mock_review = EntityReview(
        id=review_id,
        entity_type="dataset",
        entity_id=uuid4(),
        reviewer_id=TEST_USER_ID,
        status="approved",
        comments="This looks great!",
        reviewed_at=datetime.now(timezone.utc)
    )
    
    # Mock the decide_review function
    with patch('amprenta_rag.api.routers.entity_reviews.decide_review') as mock_decide_func:
        mock_decide_func.return_value = mock_review
        
        # Act
        response = client.post(
            f"/api/v1/reviews/reviews/{review_id}/decide",
            json=decision_data
        )
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(review_id)
        assert data["status"] == "approved"
        assert data["comments"] == "This looks great!"


def test_cancel_review_endpoint(client):
    """Test cancelling a review via API."""
    # Arrange
    review_id = uuid4()
    
    mock_review = EntityReview(
        id=review_id,
        entity_type="experiment",
        entity_id=uuid4(),
        reviewer_id=TEST_USER_ID,
        status="draft",
        comments="Review cancelled",
        reviewed_at=datetime.now(timezone.utc)
    )
    
    # Mock the cancel_review function
    with patch('amprenta_rag.api.routers.entity_reviews.cancel_review') as mock_cancel_func:
        mock_cancel_func.return_value = mock_review
        
        # Act
        response = client.delete(f"/api/v1/reviews/reviews/{review_id}")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["message"] == "Review cancelled successfully"


def test_create_review_endpoint_invalid_entity_type(client):
    """Test creating review with invalid entity type."""
    # Arrange
    entity_type = "invalid_type"
    entity_id = uuid4()
    
    # Act
    response = client.post(
        f"/api/v1/reviews/entities/{entity_type}/{entity_id}/reviews"
    )
    
    # Assert
    assert response.status_code == 400
    assert "Invalid entity_type" in response.json()["detail"]


def test_submit_review_endpoint_failure(client):
    """Test submitting review when service fails."""
    # Arrange
    review_id = uuid4()
    
    # Mock the submit_for_review function to return None (failure)
    with patch('amprenta_rag.api.routers.entity_reviews.submit_for_review') as mock_submit_func:
        mock_submit_func.return_value = None
        
        # Act
        response = client.post(f"/api/v1/reviews/reviews/{review_id}/submit")
        
        # Assert
        assert response.status_code == 400
        assert "Failed to submit review" in response.json()["detail"]


def test_assign_reviewer_endpoint_failure(client):
    """Test assigning reviewer when service fails."""
    # Arrange
    review_id = uuid4()
    reviewer_id = uuid4()
    
    assignment_data = {
        "reviewer_id": str(reviewer_id)
    }
    
    # Mock the assign_reviewer function to return None (failure)
    with patch('amprenta_rag.api.routers.entity_reviews.assign_reviewer') as mock_assign_func:
        mock_assign_func.return_value = None
        
        # Act
        response = client.post(
            f"/api/v1/reviews/reviews/{review_id}/assign",
            json=assignment_data
        )
        
        # Assert
        assert response.status_code == 400
        assert "Failed to assign reviewer" in response.json()["detail"]


def test_decide_review_endpoint_failure(client):
    """Test making review decision when service fails."""
    # Arrange
    review_id = uuid4()
    
    decision_data = {
        "decision": "rejected",
        "comment": "Needs more work"
    }
    
    # Mock the decide_review function to return None (failure)
    with patch('amprenta_rag.api.routers.entity_reviews.decide_review') as mock_decide_func:
        mock_decide_func.return_value = None
        
        # Act
        response = client.post(
            f"/api/v1/reviews/reviews/{review_id}/decide",
            json=decision_data
        )
        
        # Assert
        assert response.status_code == 400
        assert "Failed to make review decision" in response.json()["detail"]


def test_cancel_review_endpoint_failure(client):
    """Test cancelling review when service fails."""
    # Arrange
    review_id = uuid4()
    
    # Mock the cancel_review function to return None (failure)
    with patch('amprenta_rag.api.routers.entity_reviews.cancel_review') as mock_cancel_func:
        mock_cancel_func.return_value = None
        
        # Act
        response = client.delete(f"/api/v1/reviews/reviews/{review_id}")
        
        # Assert
        assert response.status_code == 400
        assert "Failed to cancel review" in response.json()["detail"]


def test_my_reviews_endpoint(client):
    """Test getting user's reviews via API."""
    # Arrange
    mock_reviews = [
        EntityReview(
            id=uuid4(),
            entity_type="dataset",
            entity_id=uuid4(),
            reviewer_id=TEST_USER_ID,
            status="draft",
            comments="My draft review",
            reviewed_at=datetime.now(timezone.utc)
        ),
        EntityReview(
            id=uuid4(),
            entity_type="experiment",
            entity_id=uuid4(),
            reviewer_id=TEST_USER_ID,
            status="approved",
            comments="My approved review",
            reviewed_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the get_my_reviews function
    with patch('amprenta_rag.api.routers.entity_reviews.get_my_reviews') as mock_get_func:
        mock_get_func.return_value = mock_reviews
        
        # Act
        response = client.get("/api/v1/reviews/my-reviews")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["reviews"]) == 2
        assert all(review["reviewer_id"] == str(TEST_USER_ID) for review in data["reviews"])


def test_pending_reviews_endpoint(client):
    """Test getting pending reviews via API."""
    # Arrange
    mock_reviews = [
        EntityReview(
            id=uuid4(),
            entity_type="compound",
            entity_id=uuid4(),
            reviewer_id=TEST_USER_ID,
            status="in_review",
            comments="Pending review",
            reviewed_at=datetime.now(timezone.utc)
        )
    ]
    
    # Mock the get_pending_reviews function
    with patch('amprenta_rag.api.routers.entity_reviews.get_pending_reviews') as mock_pending_func:
        mock_pending_func.return_value = mock_reviews
        
        # Act
        response = client.get("/api/v1/reviews/pending")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert len(data["reviews"]) == 1
        assert data["reviews"][0]["status"] == "in_review"


def test_status_info_endpoint(client):
    """Test getting review status information via API."""
    # Arrange
    mock_info = {
        "valid_statuses": ["draft", "submitted", "in_review", "approved", "rejected", "changes_requested"],
        "valid_transitions": {
            "draft": ["submitted"],
            "submitted": ["in_review", "draft"],
            "in_review": ["approved", "rejected", "changes_requested"]
        },
        "status_descriptions": {
            "draft": "Review is being prepared",
            "submitted": "Review submitted, waiting for reviewer assignment"
        }
    }
    
    # Mock the get_review_status_info function
    with patch('amprenta_rag.api.routers.entity_reviews.get_review_status_info') as mock_info_func:
        mock_info_func.return_value = mock_info
        
        # Act
        response = client.get("/api/v1/reviews/status-info")
        
        # Assert
        assert response.status_code == 200
        data = response.json()
        assert "valid_statuses" in data
        assert "valid_transitions" in data
        assert "status_descriptions" in data
        assert len(data["valid_statuses"]) == 6
