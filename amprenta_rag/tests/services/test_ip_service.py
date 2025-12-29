"""Tests for IP service layer."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock
from uuid import uuid4

import pytest

from amprenta_rag.services.ip_service import (
    add_patent_claim,
    create_disclosure,
    file_patent,
    get_disclosure_portfolio,
    get_entity_ip_links,
    get_patent_timeline,
    link_entity_to_disclosure,
    update_disclosure_status,
)


class TestCreateDisclosure:
    """Tests for create_disclosure."""

    def test_create_disclosure_success(self):
        """Test successful disclosure creation."""
        mock_db = MagicMock()
        
        inventors = [{"user_id": uuid4(), "contribution_percentage": 100, "is_primary": True}]
        
        result = create_disclosure(
            "Test Invention",
            "Description",
            "Chemistry",
            inventors,
            uuid4(),
            uuid4(),
            mock_db,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called

    def test_create_disclosure_with_inventors(self):
        """Test disclosure creation adds inventors."""
        mock_db = MagicMock()
        
        inventors = [
            {"user_id": uuid4(), "contribution_percentage": 60, "is_primary": True},
            {"user_id": uuid4(), "contribution_percentage": 40, "is_primary": False},
        ]
        
        result = create_disclosure(
            "Multi-Inventor",
            "Description",
            None,
            inventors,
            uuid4(),
            uuid4(),
            mock_db,
        )
        
        # Should add disclosure + 2 inventors
        assert mock_db.add.call_count >= 3


class TestUpdateDisclosureStatus:
    """Tests for update_disclosure_status."""

    def test_update_status_valid_transition(self):
        """Test valid status transition."""
        disclosure_id = uuid4()
        
        mock_disclosure = MagicMock()
        mock_disclosure.id = disclosure_id
        mock_disclosure.status = "draft"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_disclosure
        
        result = update_disclosure_status(disclosure_id, "submitted", uuid4(), mock_db)
        
        assert mock_disclosure.status == "submitted"
        assert mock_db.commit.called

    def test_update_status_invalid_transition(self):
        """Test invalid status transition raises error."""
        disclosure_id = uuid4()
        
        mock_disclosure = MagicMock()
        mock_disclosure.status = "draft"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_disclosure
        
        with pytest.raises(ValueError, match="Invalid status transition"):
            update_disclosure_status(disclosure_id, "granted", uuid4(), mock_db)


class TestFilePatent:
    """Tests for file_patent."""

    def test_file_patent_success(self):
        """Test patent filing."""
        disclosure_id = uuid4()
        
        mock_disclosure = MagicMock()
        mock_disclosure.id = disclosure_id
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_disclosure
        
        result = file_patent(
            disclosure_id,
            "US20240001234",
            "US",
            datetime(2024, 1, 1),
            mock_db,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called


class TestLinkEntity:
    """Tests for link_entity_to_disclosure."""

    def test_link_entity_success(self):
        """Test entity linking to disclosure."""
        mock_db = MagicMock()
        
        result = link_entity_to_disclosure(
            "compound",
            uuid4(),
            uuid4(),
            "embodiment",
            "Test notes",
            mock_db,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called


class TestGetDisclosurePortfolio:
    """Tests for get_disclosure_portfolio."""

    def test_get_portfolio_with_filter(self):
        """Test portfolio retrieval with status filter."""
        company_id = uuid4()
        
        mock_disclosure = MagicMock()
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.filter.return_value.order_by.return_value.all.return_value = [mock_disclosure]
        
        result = get_disclosure_portfolio(company_id, "filed", mock_db)
        
        assert isinstance(result, list)


class TestGetEntityIPLinks:
    """Tests for get_entity_ip_links."""

    def test_get_entity_links(self):
        """Test getting IP links for entity."""
        entity_id = uuid4()
        
        mock_link = MagicMock()
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_link]
        
        result = get_entity_ip_links("compound", entity_id, mock_db)
        
        assert isinstance(result, list)
        assert len(result) == 1

