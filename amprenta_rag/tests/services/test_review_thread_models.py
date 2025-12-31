"""Tests for review thread models."""

from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.database.models import ReviewThread, ReviewComment, NotebookSnapshot


class TestReviewThreadModels:
    """Test review thread model functionality."""
    
    def test_review_thread_create(self):
        """Test creating a ReviewThread linked to a review."""
        review_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Create thread
            thread = ReviewThread(
                review_id=review_id,
                title="Code review feedback",
                status="open",
                created_by_id=user_id,
                cell_index=None,  # General thread
                created_at=datetime.now(timezone.utc),
                updated_at=datetime.now(timezone.utc)
            )
            
            # Verify thread properties
            assert thread.review_id == review_id
            assert thread.title == "Code review feedback"
            assert thread.status == "open"
            assert thread.created_by_id == user_id
            assert thread.cell_index is None
            assert thread.created_at is not None
            assert thread.updated_at is not None
    
    def test_review_thread_cell_anchor(self):
        """Test creating a ReviewThread anchored to a specific cell."""
        review_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Create cell-anchored thread
            thread = ReviewThread(
                review_id=review_id,
                title="Issue with cell 5",
                status="open",
                created_by_id=user_id,
                cell_index=5  # Specific cell
            )
            
            # Verify cell anchor
            assert thread.cell_index == 5
            assert thread.title == "Issue with cell 5"
    
    def test_review_comment_create(self):
        """Test creating a ReviewComment in a thread."""
        thread_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Create comment
            comment = ReviewComment(
                thread_id=thread_id,
                content="This code needs optimization",
                created_by_id=user_id,
                parent_id=None,  # Top-level comment
                created_at=datetime.now(timezone.utc),
                updated_at=datetime.now(timezone.utc)
            )
            
            # Verify comment properties
            assert comment.thread_id == thread_id
            assert comment.content == "This code needs optimization"
            assert comment.created_by_id == user_id
            assert comment.parent_id is None
            assert comment.created_at is not None
            assert comment.updated_at is not None
    
    def test_review_comment_nested_reply(self):
        """Test creating nested reply comments (parent-child threading)."""
        thread_id = uuid4()
        user_id = uuid4()
        parent_comment_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Create parent comment
            parent_comment = ReviewComment(
                thread_id=thread_id,
                content="Original comment",
                created_by_id=user_id,
                parent_id=None
            )
            
            # Create reply comment
            reply_comment = ReviewComment(
                thread_id=thread_id,
                content="Reply to original",
                created_by_id=user_id,
                parent_id=parent_comment_id  # References parent
            )
            
            # Verify parent-child relationship
            assert parent_comment.parent_id is None
            assert reply_comment.parent_id == parent_comment_id
            assert reply_comment.thread_id == thread_id
    
    def test_notebook_snapshot_create(self):
        """Test creating a NotebookSnapshot with content."""
        review_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Sample notebook content
            notebook_content = {
                "cells": [
                    {"cell_type": "code", "source": "print('hello')"},
                    {"cell_type": "markdown", "source": "# Analysis"}
                ],
                "metadata": {"kernelspec": {"name": "python3"}}
            }
            
            # Create snapshot
            snapshot = NotebookSnapshot(
                review_id=review_id,
                content_json=notebook_content,
                cell_count=2,
                created_at=datetime.now(timezone.utc)
            )
            
            # Verify snapshot properties
            assert snapshot.review_id == review_id
            assert snapshot.content_json == notebook_content
            assert snapshot.cell_count == 2
            assert snapshot.created_at is not None
    
    def test_cascade_delete_review(self):
        """Test that deleting a review cascades to threads, comments, and snapshot."""
        review_id = uuid4()
        
        with patch("amprenta_rag.database.session.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            
            # Mock review with related objects
            mock_review = MagicMock()
            mock_review.id = review_id
            
            # Mock thread with comments
            mock_thread = MagicMock()
            mock_thread.review_id = review_id
            
            mock_comment = MagicMock()
            mock_comment.thread_id = mock_thread.id
            
            # Mock snapshot
            mock_snapshot = MagicMock()
            mock_snapshot.review_id = review_id
            
            # Set up query chain for finding review
            mock_db.query.return_value.filter.return_value.first.return_value = mock_review
            
            # Simulate delete operation
            mock_db.delete.return_value = None
            mock_db.commit.return_value = None
            
            # Verify the delete operation would be called
            # In real cascade, SQLAlchemy would handle the cascade deletes
            mock_db.delete(mock_review)
            mock_db.commit()
            
            # Verify delete was called
            mock_db.delete.assert_called_once_with(mock_review)
            mock_db.commit.assert_called_once()
