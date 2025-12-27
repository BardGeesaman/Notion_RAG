"""Unit tests for comment service utilities."""

from __future__ import annotations

from datetime import datetime, timezone
from uuid import uuid4
from unittest.mock import Mock, MagicMock

import pytest

from amprenta_rag.utils.comments import (
    update_comment,
    parse_mentions,
    resolve_mentions,
    notify_mentions,
)


def test_update_comment_success():
    """Test successful comment update."""
    comment_id = uuid4()
    user_id = uuid4()
    new_content = "Updated content"
    
    # Mock comment object
    mock_comment = Mock()
    mock_comment.id = comment_id
    mock_comment.content = "Old content"
    mock_comment.created_by_id = user_id
    
    # Mock database session
    mock_db = Mock()
    mock_query = Mock()
    mock_query.filter.return_value.first.return_value = mock_comment
    mock_db.query.return_value = mock_query
    
    result = update_comment(comment_id, new_content, user_id, mock_db)
    
    assert result == mock_comment
    assert mock_comment.content == new_content
    mock_db.commit.assert_called_once()
    mock_db.refresh.assert_called_once_with(mock_comment)


def test_update_comment_unauthorized():
    """Test update fails when user is not the author."""
    comment_id = uuid4()
    user_id = uuid4()
    other_user_id = uuid4()
    
    # Mock database session
    mock_db = Mock()
    mock_query = Mock()
    mock_query.filter.return_value.first.return_value = None  # Not found for this user
    mock_db.query.return_value = mock_query
    
    result = update_comment(comment_id, "New content", user_id, mock_db)
    
    assert result is None
    mock_db.commit.assert_not_called()


def test_parse_mentions_single():
    """Test parsing a single @mention."""
    content = "Hello @john, how are you?"
    mentions = parse_mentions(content)
    
    assert mentions == ["john"]


def test_parse_mentions_multiple():
    """Test parsing multiple @mentions."""
    content = "Hey @alice and @bob, meet @charlie_123"
    mentions = parse_mentions(content)
    
    assert set(mentions) == {"alice", "bob", "charlie_123"}
    assert len(mentions) == 3


def test_parse_mentions_empty():
    """Test parsing content with no mentions."""
    content = "This is a comment without any mentions."
    mentions = parse_mentions(content)
    
    assert mentions == []


def test_resolve_mentions():
    """Test resolving usernames to user IDs."""
    usernames = ["alice", "bob"]
    user_id_1 = uuid4()
    user_id_2 = uuid4()
    
    # Mock User objects
    mock_user_1 = Mock()
    mock_user_1.id = user_id_1
    mock_user_2 = Mock()
    mock_user_2.id = user_id_2
    
    # Mock database session
    mock_db = Mock()
    mock_query = Mock()
    mock_query.filter.return_value.all.return_value = [mock_user_1, mock_user_2]
    mock_db.query.return_value = mock_query
    
    result = resolve_mentions(usernames, mock_db)
    
    assert result == [user_id_1, user_id_2]


def test_notify_mentions():
    """Test notifying mentioned users."""
    comment_id = uuid4()
    author_id = uuid4()
    mentioned_user_1 = uuid4()
    mentioned_user_2 = uuid4()
    
    # Mock comment
    mock_comment = Mock()
    mock_comment.id = comment_id
    mock_comment.created_by_id = author_id
    mock_comment.entity_type = "dataset"
    mock_comment.entity_id = uuid4()
    mock_comment.content = "Test comment content"
    
    mock_db = Mock()
    
    # Mock log_activity to avoid actual activity logging
    from unittest.mock import patch
    with patch("amprenta_rag.utils.comments.log_activity") as mock_log:
        notify_mentions(mock_comment, [mentioned_user_1, mentioned_user_2], mock_db)
        
        # Should be called twice (once per mentioned user)
        assert mock_log.call_count == 2


def test_notify_mentions_excludes_author():
    """Test that notification excludes the comment author."""
    comment_id = uuid4()
    author_id = uuid4()
    
    # Mock comment
    mock_comment = Mock()
    mock_comment.id = comment_id
    mock_comment.created_by_id = author_id
    mock_comment.entity_type = "experiment"
    mock_comment.entity_id = uuid4()
    mock_comment.content = "Self mention @me"
    
    mock_db = Mock()
    
    # Mock log_activity
    from unittest.mock import patch
    with patch("amprenta_rag.utils.comments.log_activity") as mock_log:
        # Include author in mentioned users
        notify_mentions(mock_comment, [author_id], mock_db)
        
        # Should not be called (author excluded)
        mock_log.assert_not_called()
