"""Tests for share link service."""

from __future__ import annotations

from datetime import datetime, timedelta
from unittest.mock import MagicMock
from uuid import uuid4

import pytest

from amprenta_rag.services.share_link_service import (
    cleanup_expired_links,
    generate_share_link,
    get_share_link_stats,
    get_user_share_links,
    revoke_share_link,
    validate_share_link,
)


class TestGenerateShareLink:
    """Tests for generate_share_link."""

    def test_generate_success(self):
        """Test successful share link generation."""
        mock_db = MagicMock()
        
        result = generate_share_link(
            "/voila/dashboard.ipynb",
            {"experiment_id": "123"},
            uuid4(),
            uuid4(),
            24,
            10,
            "view",
            mock_db,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called

    def test_token_is_secure(self):
        """Test that generated token is 64 chars."""
        mock_db = MagicMock()
        
        # Mock the ShareLink to capture token
        added_links = []
        mock_db.add.side_effect = lambda obj: added_links.append(obj)
        
        generate_share_link(
            "/voila/test.ipynb",
            None,
            uuid4(),
            uuid4(),
            1,
            None,
            "view",
            mock_db,
        )
        
        # Token should be 64 URL-safe chars from token_urlsafe(48)
        assert len(added_links) > 0
        token = added_links[0].token
        assert len(token) == 64


class TestValidateShareLink:
    """Tests for validate_share_link."""

    def test_validate_success(self):
        """Test successful validation."""
        mock_link = MagicMock()
        mock_link.is_active = True
        mock_link.expires_at = datetime.utcnow() + timedelta(hours=1)
        mock_link.max_views = None
        mock_link.view_count = 0
        mock_link.dashboard_path = "/voila/dashboard.ipynb"
        mock_link.context_json = '{"key": "value"}'
        mock_link.permissions = "view"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_link
        
        result = validate_share_link("test_token", mock_db)
        
        assert result is not None
        assert result["dashboard_path"] == "/voila/dashboard.ipynb"
        assert mock_db.commit.called  # View count incremented

    def test_validate_expired(self):
        """Test validation of expired link."""
        mock_link = MagicMock()
        mock_link.is_active = True
        mock_link.expires_at = datetime.utcnow() - timedelta(hours=1)  # Expired
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_link
        
        result = validate_share_link("expired_token", mock_db)
        
        assert result is None
        assert mock_link.is_active is False  # Should be deactivated

    def test_validate_max_views_exceeded(self):
        """Test validation with max views exceeded."""
        mock_link = MagicMock()
        mock_link.is_active = True
        mock_link.expires_at = datetime.utcnow() + timedelta(hours=1)
        mock_link.max_views = 5
        mock_link.view_count = 5  # At limit
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_link
        
        result = validate_share_link("maxed_token", mock_db)
        
        assert result is None


class TestRevokeShareLink:
    """Tests for revoke_share_link."""

    def test_revoke_success(self):
        """Test successful revocation."""
        link_id = uuid4()
        user_id = uuid4()
        
        mock_link = MagicMock()
        mock_link.id = link_id
        mock_link.created_by_id = user_id
        mock_link.is_active = True
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_link
        
        result = revoke_share_link(link_id, user_id, mock_db)
        
        assert result is True
        assert mock_link.is_active is False


class TestCleanupExpiredLinks:
    """Tests for cleanup_expired_links."""

    def test_cleanup(self):
        """Test cleanup of expired links."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.update.return_value = 3
        
        count = cleanup_expired_links(mock_db)
        
        assert count == 3
        assert mock_db.commit.called

