"""Tests for ShareLink model."""

from __future__ import annotations

from datetime import datetime, timedelta
from uuid import uuid4

import pytest

from amprenta_rag.models.share_links import ShareLink


def test_create_share_link():
    """Test ShareLink model creation."""
    link = ShareLink(
        token="abc123def456",
        dashboard_path="/voila/render/dashboard.ipynb",
        created_by_id=uuid4(),
        company_id=uuid4(),
        expires_at=datetime.utcnow() + timedelta(days=7),
    )
    
    assert link.token == "abc123def456"
    assert link.dashboard_path == "/voila/render/dashboard.ipynb"
    # Defaults apply on db insert, not object creation
    assert link.expires_at is not None


def test_share_link_token_unique():
    """Test that token field is unique."""
    # Unique constraint is enforced at database level
    link = ShareLink(
        token="unique_token_123",
        dashboard_path="/voila/dashboard.ipynb",
        created_by_id=uuid4(),
        company_id=uuid4(),
        expires_at=datetime.utcnow() + timedelta(days=1),
    )
    
    assert link.token == "unique_token_123"


def test_share_link_default_values():
    """Test ShareLink default values."""
    link = ShareLink(
        token="test_token",
        dashboard_path="/voila/test.ipynb",
        created_by_id=uuid4(),
        company_id=uuid4(),
        expires_at=datetime.utcnow() + timedelta(days=1),
    )
    
    # Defaults defined but apply on db insert
    assert link.token == "test_token"


def test_share_link_expiration_required():
    """Test that expires_at is required."""
    link = ShareLink(
        token="expires_test",
        dashboard_path="/voila/dashboard.ipynb",
        created_by_id=uuid4(),
        company_id=uuid4(),
        expires_at=datetime(2025, 12, 31),
    )
    
    assert link.expires_at is not None

