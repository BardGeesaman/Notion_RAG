"""API tests for notebook review threads and diff functionality."""

from datetime import datetime, timezone
from unittest.mock import MagicMock, patch, ANY
from uuid import uuid4

import pytest
from fastapi import status
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.models import User, NotebookReview, ReviewThread, ReviewComment
from amprenta_rag.api.dependencies import get_current_user, get_database_session


class TestReviewThreadsAPI:
    """Test cases for review threads API endpoints."""

    @pytest.fixture
    def mock_user(self):
        """Mock authenticated user."""
        return User(
            id=uuid4(),
            email="test@example.com",
            username="testuser",
            role="researcher",
        )

    @pytest.fixture
    def mock_review(self):
        """Mock notebook review."""
        return NotebookReview(
            id=uuid4(),
            notebook_path="path/to/test_notebook.ipynb",
            version_hash="abc123def456",
            status="pending",
            created_at=datetime.now(timezone.utc),
        )

    @pytest.fixture
    def mock_thread(self, mock_review, mock_user):
        """Mock review thread."""
        return ReviewThread(
            id=uuid4(),
            review_id=mock_review.id,
            title="Initial feedback on introduction",
            status="open",
            created_by_id=mock_user.id,
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

    @pytest.fixture
    def mock_comment(self, mock_thread, mock_user):
        """Mock thread comment."""
        return ReviewComment(
            id=uuid4(),
            thread_id=mock_thread.id,
            content="This section needs more detail.",
            created_by_id=mock_user.id,
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

    @pytest.fixture
    def client(self, mock_user):
        """Test client with mocked authentication."""
        def mock_get_user():
            return mock_user

        def mock_get_db():
            return MagicMock()

        app.dependency_overrides[get_current_user] = mock_get_user
        app.dependency_overrides[get_database_session] = mock_get_db

        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()

    def test_create_thread_success(self, client, mock_review, mock_thread, mock_user):
        """Test successful creation of a review thread."""
        with patch("amprenta_rag.api.routers.review_threads.create_thread") as mock_create:
            mock_create.return_value = mock_thread

            response = client.post(
                f"/api/v1/reviews/{mock_review.id}/threads",
                json={
                    "title": "Initial feedback on introduction",
                    "cell_index": None,
                },
            )

            assert response.status_code == status.HTTP_201_CREATED
            data = response.json()
            assert data["id"] == str(mock_thread.id)
            assert data["title"] == "Initial feedback on introduction"
            assert data["status"] == "open"
            assert data["created_by_id"] == str(mock_user.id)
            assert data["comments"] == []

            mock_create.assert_called_once_with(
                db=ANY,
                review_id=mock_review.id,
                title="Initial feedback on introduction",
                created_by_id=mock_user.id,
                cell_index=None,
            )

    def test_create_thread_with_cell_index(self, client, mock_review, mock_user):
        """Test creating a thread anchored to a specific cell."""
        mock_thread = ReviewThread(
            id=uuid4(),
            review_id=mock_review.id,
            title="Cell 3 needs clarification",
            cell_index=3,
            status="open",
            created_by_id=mock_user.id,
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

        with patch("amprenta_rag.api.routers.review_threads.create_thread") as mock_create:
            mock_create.return_value = mock_thread

            response = client.post(
                f"/api/v1/reviews/{mock_review.id}/threads",
                json={
                    "title": "Cell 3 needs clarification",
                    "cell_index": 3,
                },
            )

            assert response.status_code == status.HTTP_201_CREATED
            data = response.json()
            assert data["cell_index"] == 3
            assert data["title"] == "Cell 3 needs clarification"

            mock_create.assert_called_once_with(
                db=ANY,
                review_id=mock_review.id,
                title="Cell 3 needs clarification",
                created_by_id=mock_user.id,
                cell_index=3,
            )

    def test_list_threads_empty(self, client, mock_review):
        """Test listing threads when none exist."""
        with patch("amprenta_rag.api.routers.review_threads.get_threads") as mock_get:
            mock_get.return_value = []

            response = client.get(f"/api/v1/reviews/{mock_review.id}/threads")

            assert response.status_code == status.HTTP_200_OK
            data = response.json()
            assert data == []

            mock_get.assert_called_once_with(
                db=ANY,
                review_id=mock_review.id,
            )

    def test_list_threads_with_comments(self, client, mock_review, mock_thread, mock_comment, mock_user):
        """Test listing threads with nested comments."""
        # Set up thread with comments
        mock_thread.comments = [mock_comment]

        with patch("amprenta_rag.api.routers.review_threads.get_threads") as mock_get:
            mock_get.return_value = [mock_thread]

            response = client.get(f"/api/v1/reviews/{mock_review.id}/threads")

            assert response.status_code == status.HTTP_200_OK
            data = response.json()
            assert len(data) == 1

            thread_data = data[0]
            assert thread_data["id"] == str(mock_thread.id)
            assert thread_data["title"] == mock_thread.title
            assert len(thread_data["comments"]) == 1

            comment_data = thread_data["comments"][0]
            assert comment_data["id"] == str(mock_comment.id)
            assert comment_data["content"] == mock_comment.content
            assert comment_data["created_by_id"] == str(mock_user.id)

    def test_add_comment_success(self, client, mock_thread, mock_comment, mock_user):
        """Test adding a comment to a thread."""
        with patch("amprenta_rag.api.routers.review_threads.add_comment") as mock_add:
            mock_add.return_value = mock_comment

            response = client.post(
                f"/api/v1/threads/{mock_thread.id}/comments",
                json={
                    "content": "This section needs more detail.",
                    "parent_id": None,
                },
            )

            assert response.status_code == status.HTTP_201_CREATED
            data = response.json()
            assert data["id"] == str(mock_comment.id)
            assert data["content"] == "This section needs more detail."
            assert data["thread_id"] == str(mock_thread.id)
            assert data["parent_id"] is None
            assert data["created_by_id"] == str(mock_user.id)

            mock_add.assert_called_once_with(
                db=ANY,
                thread_id=mock_thread.id,
                content="This section needs more detail.",
                created_by_id=mock_user.id,
                parent_id=None,
            )

    def test_add_nested_reply(self, client, mock_thread, mock_user):
        """Test adding a nested reply to a comment."""
        parent_comment_id = uuid4()
        reply_comment = ReviewComment(
            id=uuid4(),
            thread_id=mock_thread.id,
            parent_id=parent_comment_id,
            content="Good point, I'll revise that section.",
            created_by_id=mock_user.id,
            created_at=datetime.now(timezone.utc),
            updated_at=datetime.now(timezone.utc),
        )

        with patch("amprenta_rag.api.routers.review_threads.add_comment") as mock_add:
            mock_add.return_value = reply_comment

            response = client.post(
                f"/api/v1/threads/{mock_thread.id}/comments",
                json={
                    "content": "Good point, I'll revise that section.",
                    "parent_id": str(parent_comment_id),
                },
            )

            assert response.status_code == status.HTTP_201_CREATED
            data = response.json()
            assert data["parent_id"] == str(parent_comment_id)
            assert data["content"] == "Good point, I'll revise that section."

            mock_add.assert_called_once_with(
                db=ANY,
                thread_id=mock_thread.id,
                content="Good point, I'll revise that section.",
                created_by_id=mock_user.id,
                parent_id=parent_comment_id,
            )

    def test_update_thread_status(self, client, mock_thread, mock_user):
        """Test updating thread status to resolved."""
        resolved_thread = ReviewThread(
            id=mock_thread.id,
            review_id=mock_thread.review_id,
            title=mock_thread.title,
            status="resolved",
            created_by_id=mock_user.id,
            created_at=mock_thread.created_at,
            updated_at=datetime.now(timezone.utc),
        )

        with patch("amprenta_rag.api.routers.review_threads.resolve_thread") as mock_resolve:
            mock_resolve.return_value = resolved_thread

            response = client.patch(
                f"/api/v1/threads/{mock_thread.id}",
                json={"status": "resolved"},
            )

            assert response.status_code == status.HTTP_200_OK
            data = response.json()
            assert data["id"] == str(mock_thread.id)
            assert data["status"] == "resolved"

            mock_resolve.assert_called_once_with(
                db=ANY,
                thread_id=mock_thread.id,
                status="resolved",
            )

    def test_get_diff_success(self, client, mock_review, mock_user):
        """Test getting notebook diff between snapshot and current version."""
        mock_diff_result = {
            "added": [
                {"index": 0, "cell_type": "code", "source_preview": "print('new cell')"}
            ],
            "removed": [
                {"index": 1, "cell_type": "markdown", "source_preview": "# Old section"}
            ],
            "modified": [
                {
                    "index": 2,
                    "cell_type": "code",
                    "old_source_preview": "x = 1",
                    "new_source_preview": "x = 2",
                }
            ],
            "unchanged": [
                {"index": 3, "cell_type": "markdown", "source_preview": "# Conclusion"}
            ],
        }

        with patch("amprenta_rag.api.routers.review_threads.get_review_diff") as mock_diff:
            mock_diff.return_value = mock_diff_result

            response = client.get(f"/api/v1/reviews/{mock_review.id}/diff")

            assert response.status_code == status.HTTP_200_OK
            data = response.json()

            assert len(data["added"]) == 1
            assert data["added"][0]["source_preview"] == "print('new cell')"

            assert len(data["removed"]) == 1
            assert data["removed"][0]["source_preview"] == "# Old section"

            assert len(data["modified"]) == 1
            assert data["modified"][0]["old_source_preview"] == "x = 1"
            assert data["modified"][0]["new_source_preview"] == "x = 2"

            assert len(data["unchanged"]) == 1
            assert data["unchanged"][0]["source_preview"] == "# Conclusion"

            # Just verify the function was called with the right review ID
            # Don't check exact parameters since the notebook_path comes from DB
            assert mock_diff.called
            call_args = mock_diff.call_args[1]  # Get keyword arguments
            assert call_args["review_id"] == mock_review.id

    def test_create_thread_invalid_review(self, client):
        """Test creating thread with non-existent review ID."""
        with patch("amprenta_rag.api.routers.review_threads.create_thread") as mock_create:
            mock_create.side_effect = ValueError("NotebookReview with ID not found")

            response = client.post(
                f"/api/v1/reviews/{uuid4()}/threads",
                json={"title": "Test thread"},
            )

            assert response.status_code == status.HTTP_404_NOT_FOUND
            assert "not found" in response.json()["detail"]

    def test_add_comment_invalid_thread(self, client):
        """Test adding comment to non-existent thread ID."""
        with patch("amprenta_rag.api.routers.review_threads.add_comment") as mock_add:
            mock_add.side_effect = ValueError("ReviewThread with ID not found")

            response = client.post(
                f"/api/v1/threads/{uuid4()}/comments",
                json={"content": "Test comment"},
            )

            assert response.status_code == status.HTTP_404_NOT_FOUND
            assert "not found" in response.json()["detail"]

    def test_update_thread_invalid_status(self, client, mock_thread):
        """Test updating thread with invalid status."""
        response = client.patch(
            f"/api/v1/threads/{mock_thread.id}",
            json={"status": "invalid_status"},
        )

        assert response.status_code == status.HTTP_400_BAD_REQUEST
        assert "Status must be" in response.json()["detail"]

    def test_get_diff_review_not_found(self, client):
        """Test getting diff for non-existent review."""
        mock_db = client.app.dependency_overrides[get_database_session]()
        mock_db.query.return_value.filter.return_value.first.return_value = None

        response = client.get(f"/api/v1/reviews/{uuid4()}/diff")

        assert response.status_code == status.HTTP_404_NOT_FOUND
        assert "not found" in response.json()["detail"]
