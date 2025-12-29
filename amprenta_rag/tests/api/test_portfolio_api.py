"""Tests for portfolio API endpoints."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestPortfolioSummary:
    """Tests for GET /api/v1/portfolio/summary endpoint."""

    @patch("amprenta_rag.api.routers.portfolio.get_portfolio_summary")
    def test_get_summary_success(self, mock_get_summary):
        """Test successful portfolio summary."""
        mock_get_summary.return_value = {
            "total_compounds": 150,
            "scaffold_count": 25,
            "date_from": "2024-01-01T00:00:00",
            "date_to": "2024-12-31T23:59:59",
        }
        
        response = client.get("/api/v1/portfolio/summary")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total_compounds"] == 150
        assert data["scaffold_count"] == 25


class TestADMETSummary:
    """Tests for GET /api/v1/portfolio/admet endpoint."""

    @patch("amprenta_rag.api.routers.portfolio.get_admet_summary")
    def test_get_admet_success(self, mock_get_admet):
        """Test successful ADMET summary."""
        mock_get_admet.return_value = {
            "green": 80,
            "yellow": 50,
            "red": 20,
            "unknown": 0,
        }
        
        response = client.get("/api/v1/portfolio/admet")
        
        assert response.status_code == 200
        data = response.json()
        assert data["green"] == 80
        assert data["yellow"] == 50
        assert data["red"] == 20


class TestSARGaps:
    """Tests for GET /api/v1/portfolio/gaps endpoint."""

    @patch("amprenta_rag.api.routers.portfolio.get_sar_gaps")
    def test_get_gaps_success(self, mock_get_gaps):
        """Test successful SAR gaps retrieval."""
        mock_get_gaps.return_value = [
            {"scaffold_id": "scf_001", "compound_count": 2, "example_smiles": "c1ccccc1"},
        ]
        
        response = client.get("/api/v1/portfolio/gaps")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["scaffold_id"] == "scf_001"


class TestRecommendations:
    """Tests for GET /api/v1/portfolio/recommendations endpoint."""

    @patch("amprenta_rag.api.routers.portfolio.get_recommendations")
    def test_get_recommendations_success(self, mock_get_recs):
        """Test successful recommendations retrieval."""
        mock_get_recs.return_value = [
            {"compound_id": "CMPD001", "smiles": "CCO", "score": 0.95, "reason": "High score"},
        ]
        
        response = client.get("/api/v1/portfolio/recommendations?limit=5")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["compound_id"] == "CMPD001"

    @patch("amprenta_rag.api.routers.portfolio.get_recommendations")
    def test_get_recommendations_empty(self, mock_get_recs):
        """Test recommendations when empty."""
        mock_get_recs.return_value = []
        
        response = client.get("/api/v1/portfolio/recommendations")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []

