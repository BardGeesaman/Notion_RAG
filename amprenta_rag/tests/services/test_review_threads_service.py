"""Tests for review threads and notebook diff services."""

import json
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch, mock_open
from uuid import uuid4

from amprenta_rag.services.review_threads import (
    create_thread,
    get_threads,
    add_comment,
    resolve_thread,
    reopen_thread,
    delete_thread
)
from amprenta_rag.services.notebook_diff import (
    capture_snapshot,
    compute_cell_diff
)


class TestReviewThreadsService:
    """Test review threads service functionality."""
    
    def test_create_thread_success(self):
        """Test creating a review thread successfully."""
        review_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.services.review_threads.ReviewThread") as mock_thread_class:
            mock_db = MagicMock()
            mock_thread = MagicMock()
            mock_thread.id = uuid4()
            mock_thread_class.return_value = mock_thread
            
            result = create_thread(
                db=mock_db,
                review_id=review_id,
                title="Code review feedback",
                created_by_id=user_id,
                cell_index=None
            )
            
            # Verify thread creation
            mock_thread_class.assert_called_once_with(
                review_id=review_id,
                title="Code review feedback",
                created_by_id=user_id,
                cell_index=None,
                status="open"
            )
            mock_db.add.assert_called_once_with(mock_thread)
            mock_db.commit.assert_called_once()
            mock_db.refresh.assert_called_once_with(mock_thread)
            assert result == mock_thread
    
    def test_create_thread_with_cell_anchor(self):
        """Test creating a thread anchored to a specific cell."""
        review_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.services.review_threads.ReviewThread") as mock_thread_class:
            mock_db = MagicMock()
            mock_thread = MagicMock()
            mock_thread_class.return_value = mock_thread
            
            create_thread(
                db=mock_db,
                review_id=review_id,
                title="Issue with cell 5",
                created_by_id=user_id,
                cell_index=5
            )
            
            # Verify cell_index is passed correctly
            mock_thread_class.assert_called_once_with(
                review_id=review_id,
                title="Issue with cell 5",
                created_by_id=user_id,
                cell_index=5,
                status="open"
            )
    
    def test_get_threads_empty(self):
        """Test getting threads when none exist."""
        review_id = uuid4()
        
        mock_db = MagicMock()
        mock_db.query.return_value.options.return_value.filter.return_value.order_by.return_value.all.return_value = []
        
        result = get_threads(mock_db, review_id)
        
        assert result == []
        mock_db.query.assert_called_once()
    
    def test_get_threads_with_comments(self):
        """Test getting threads with comments eager loaded."""
        review_id = uuid4()
        
        mock_thread1 = MagicMock()
        mock_thread1.id = uuid4()
        mock_thread2 = MagicMock()
        mock_thread2.id = uuid4()
        
        mock_db = MagicMock()
        mock_db.query.return_value.options.return_value.filter.return_value.order_by.return_value.all.return_value = [
            mock_thread1, mock_thread2
        ]
        
        result = get_threads(mock_db, review_id)
        
        assert len(result) == 2
        assert result[0] == mock_thread1
        assert result[1] == mock_thread2
    
    def test_add_comment_to_thread(self):
        """Test adding a comment to a thread."""
        thread_id = uuid4()
        user_id = uuid4()
        
        with patch("amprenta_rag.services.review_threads.ReviewComment") as mock_comment_class:
            mock_db = MagicMock()
            mock_comment = MagicMock()
            mock_comment.id = uuid4()
            mock_comment_class.return_value = mock_comment
            
            result = add_comment(
                db=mock_db,
                thread_id=thread_id,
                content="This needs optimization",
                created_by_id=user_id,
                parent_id=None
            )
            
            # Verify comment creation
            mock_comment_class.assert_called_once_with(
                thread_id=thread_id,
                content="This needs optimization",
                created_by_id=user_id,
                parent_id=None
            )
            mock_db.add.assert_called_once_with(mock_comment)
            mock_db.commit.assert_called_once()
            mock_db.refresh.assert_called_once_with(mock_comment)
            assert result == mock_comment
    
    def test_add_nested_reply(self):
        """Test adding a nested reply to a comment."""
        thread_id = uuid4()
        user_id = uuid4()
        parent_id = uuid4()
        
        with patch("amprenta_rag.services.review_threads.ReviewComment") as mock_comment_class:
            mock_db = MagicMock()
            mock_comment = MagicMock()
            mock_comment_class.return_value = mock_comment
            
            add_comment(
                db=mock_db,
                thread_id=thread_id,
                content="I agree with the above",
                created_by_id=user_id,
                parent_id=parent_id
            )
            
            # Verify parent_id is passed correctly
            mock_comment_class.assert_called_once_with(
                thread_id=thread_id,
                content="I agree with the above",
                created_by_id=user_id,
                parent_id=parent_id
            )
    
    def test_resolve_thread(self):
        """Test resolving a thread."""
        thread_id = uuid4()
        
        mock_thread = MagicMock()
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_thread
        
        result = resolve_thread(mock_db, thread_id, "resolved")
        
        assert mock_thread.status == "resolved"
        mock_db.commit.assert_called_once()
        mock_db.refresh.assert_called_once_with(mock_thread)
        assert result == mock_thread
    
    def test_reopen_thread(self):
        """Test reopening a resolved thread."""
        thread_id = uuid4()
        
        mock_thread = MagicMock()
        mock_thread.status = "resolved"
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_thread
        
        result = reopen_thread(mock_db, thread_id)
        
        assert mock_thread.status == "open"
        mock_db.commit.assert_called_once()
        mock_db.refresh.assert_called_once_with(mock_thread)
        assert result == mock_thread


class TestNotebookDiffService:
    """Test notebook diff service functionality."""
    
    @patch("amprenta_rag.services.notebook_diff.load_registry")
    @patch("amprenta_rag.services.notebook_diff.resolve_notebook_path")
    @patch("amprenta_rag.services.notebook_diff.nbformat")
    @patch("builtins.open", new_callable=mock_open)
    def test_capture_snapshot_success(self, mock_file, mock_nbformat, mock_resolve_path, mock_load_registry):
        """Test capturing a notebook snapshot successfully."""
        review_id = uuid4()
        notebook_path = "analysis/example.ipynb"
        
        # Mock registry and path resolution
        mock_template = {"notebook_path": notebook_path}
        mock_load_registry.return_value = [mock_template]
        mock_resolve_path.return_value = MagicMock(exists=lambda: True)
        
        # Mock notebook content
        mock_notebook = MagicMock()
        mock_nbformat.read.return_value = mock_notebook
        mock_nbformat.v4.writes.return_value = json.dumps({
            "cells": [
                {"cell_type": "code", "source": "print('hello')"},
                {"cell_type": "markdown", "source": "# Analysis"}
            ],
            "metadata": {}
        })
        
        # Mock database
        with patch("amprenta_rag.services.notebook_diff.NotebookSnapshot") as mock_snapshot_class:
            mock_db = MagicMock()
            mock_snapshot = MagicMock()
            mock_snapshot.id = uuid4()
            mock_snapshot_class.return_value = mock_snapshot
            
            result = capture_snapshot(mock_db, review_id, notebook_path)
            
            # Verify snapshot creation
            mock_snapshot_class.assert_called_once()
            call_args = mock_snapshot_class.call_args[1]
            assert call_args["review_id"] == review_id
            assert call_args["cell_count"] == 2
            assert "cells" in call_args["content_json"]
            
            mock_db.add.assert_called_once_with(mock_snapshot)
            mock_db.commit.assert_called_once()
            mock_db.refresh.assert_called_once_with(mock_snapshot)
            assert result == mock_snapshot
    
    def test_compute_cell_diff_modifications(self):
        """Test computing differences between notebook cells."""
        old_cells = [
            {"cell_type": "code", "source": "print('hello')"},
            {"cell_type": "code", "source": "x = 1"}
        ]
        
        new_cells = [
            {"cell_type": "code", "source": "print('hello')"},  # Unchanged
            {"cell_type": "code", "source": "x = 1"},           # Unchanged
            {"cell_type": "code", "source": "y = 2"}            # Added
        ]
        
        result = compute_cell_diff(old_cells, new_cells)
        
        # Verify diff structure
        assert "added" in result
        assert "removed" in result
        assert "modified" in result
        assert "unchanged" in result
        
        # Verify we have changes detected
        assert len(result["added"]) >= 1, f"Expected at least 1 added cell, got {len(result['added'])}"
        assert len(result["unchanged"]) >= 1, f"Expected at least 1 unchanged cell, got {len(result['unchanged'])}"
        
        # Verify we can find the new cell
        added_cells = result["added"]
        new_cell_found = any("y = 2" in cell["source_preview"] for cell in added_cells)
        assert new_cell_found, "Should find the newly added cell with 'y = 2'"
        
        # Verify we can find an unchanged cell
        unchanged_cells = result["unchanged"]
        unchanged_cell_found = any("x = 1" in cell["source_preview"] for cell in unchanged_cells)
        assert unchanged_cell_found, "Should find the unchanged cell with 'x = 1'"
