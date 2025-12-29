"""Tests for portfolio service."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock
from uuid import uuid4

import pytest

from amprenta_rag.services.portfolio_service import (
    get_portfolio_summary,
    get_admet_summary,
    get_sar_gaps,
    get_recommendations,
)


class TestPortfolioSummary:
    """Tests for get_portfolio_summary."""

    def test_get_summary_with_compounds(self):
        """Test summary with compounds in database."""
        mock_db = MagicMock()
        mock_db.query.return_value.count.return_value = 150
        
        mock_first = MagicMock()
        mock_first.created_at = datetime(2024, 1, 1)
        
        mock_last = MagicMock()
        mock_last.created_at = datetime(2024, 12, 31)
        
        mock_db.query.return_value.order_by.return_value.first.side_effect = [mock_first, mock_last]
        
        result = get_portfolio_summary(mock_db)
        
        assert result["total_compounds"] == 150
        assert result["date_from"] is not None
        assert result["date_to"] is not None

    def test_get_summary_empty_database(self):
        """Test summary with no compounds."""
        mock_db = MagicMock()
        mock_db.query.return_value.count.return_value = 0
        mock_db.query.return_value.order_by.return_value.first.return_value = None
        
        result = get_portfolio_summary(mock_db)
        
        assert result["total_compounds"] == 0


class TestADMETSummary:
    """Tests for get_admet_summary."""

    def test_get_admet_summary_placeholder(self):
        """Test ADMET summary returns structure."""
        mock_db = MagicMock()
        mock_db.query.return_value.count.return_value = 100
        
        result = get_admet_summary(mock_db)
        
        assert "green" in result
        assert "yellow" in result
        assert "red" in result
        assert "unknown" in result


class TestSARGaps:
    """Tests for get_sar_gaps."""

    def test_get_sar_gaps_returns_list(self):
        """Test SAR gaps returns list."""
        mock_db = MagicMock()
        
        result = get_sar_gaps(mock_db, min_compounds=3)
        
        assert isinstance(result, list)


class TestRecommendations:
    """Tests for get_recommendations."""

    def test_get_recommendations_with_compounds(self):
        """Test recommendations with compounds."""
        mock_comp = MagicMock()
        mock_comp.compound_id = "CMPD001"
        mock_comp.id = uuid4()
        mock_comp.canonical_smiles = "CCO"
        mock_comp.smiles = "CCO"
        
        mock_db = MagicMock()
        mock_db.query.return_value.order_by.return_value.limit.return_value.all.return_value = [mock_comp]
        
        result = get_recommendations(mock_db, limit=5)
        
        assert len(result) == 1
        assert result[0]["compound_id"] == "CMPD001"
        assert result[0]["smiles"] == "CCO"

    def test_get_recommendations_empty(self):
        """Test recommendations with no compounds."""
        mock_db = MagicMock()
        mock_db.query.return_value.order_by.return_value.limit.return_value.all.return_value = []
        
        result = get_recommendations(mock_db, limit=5)
        
        assert result == []

