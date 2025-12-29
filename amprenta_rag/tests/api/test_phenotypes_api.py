"""
Unit tests for phenotypes API endpoints.

Tests HPO phenotype mapping with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestGenesForHPO:
    """Tests for GET /api/phenotypes/{hpo_id}/genes endpoint."""

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_genes_for_hpo_success(self, mock_get_mapper):
        """Test successful gene retrieval for HPO term."""
        hpo_id = "HP:0000001"
        
        # Mock mapper
        mock_mapper = MagicMock()
        mock_mapper.get_genes_for_hpo.return_value = ["SOD1", "TARDBP", "FUS"]
        mock_get_mapper.return_value = mock_mapper
        
        response = client.get(f"/api/phenotypes/{hpo_id}/genes")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 3
        assert "SOD1" in data
        assert "TARDBP" in data
        assert "FUS" in data

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_genes_for_hpo_empty(self, mock_get_mapper):
        """Test gene retrieval for HPO term with no associations."""
        hpo_id = "HP:9999999"
        
        # Mock mapper returning empty list
        mock_mapper = MagicMock()
        mock_mapper.get_genes_for_hpo.return_value = []
        mock_get_mapper.return_value = mock_mapper
        
        response = client.get(f"/api/phenotypes/{hpo_id}/genes")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 0

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_genes_for_hpo_error(self, mock_get_mapper):
        """Test gene retrieval with lookup error."""
        hpo_id = "HP:0000001"
        
        # Mock mapper raising exception
        mock_mapper = MagicMock()
        mock_mapper.get_genes_for_hpo.side_effect = RuntimeError("Lookup failed")
        mock_get_mapper.return_value = mock_mapper
        
        response = client.get(f"/api/phenotypes/{hpo_id}/genes")
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()


class TestExpandQuery:
    """Tests for POST /api/phenotypes/expand-query endpoint."""

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_expand_query_success(self, mock_get_mapper):
        """Test successful query expansion with HPO IDs."""
        # Mock mapper
        mock_mapper = MagicMock()
        mock_mapper.expand_query.return_value = {
            "hpo_ids": ["HP:0000001", "HP:0000002"],
            "genes": ["SOD1", "TARDBP", "FUS", "C9orf72"],
            "gene_count": 4,
        }
        mock_get_mapper.return_value = mock_mapper
        
        response = client.post(
            "/api/phenotypes/expand-query",
            json={"query": "Patient has HP:0000002 and HP:0000001"},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["hpo_ids"] == ["HP:0000001", "HP:0000002"]
        assert isinstance(data["genes"], list)
        assert data["gene_count"] == 4
        assert "SOD1" in data["genes"]

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_expand_query_no_hpo_ids(self, mock_get_mapper):
        """Test query expansion with no HPO IDs found."""
        # Mock mapper
        mock_mapper = MagicMock()
        mock_mapper.expand_query.return_value = {
            "hpo_ids": [],
            "genes": [],
            "gene_count": 0,
        }
        mock_get_mapper.return_value = mock_mapper
        
        response = client.post(
            "/api/phenotypes/expand-query",
            json={"query": "Patient has symptoms"},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["hpo_ids"] == []
        assert data["genes"] == []
        assert data["gene_count"] == 0

    @patch("amprenta_rag.api.routers.phenotypes.get_mapper")
    def test_expand_query_error(self, mock_get_mapper):
        """Test query expansion with error."""
        # Mock mapper raising exception
        mock_mapper = MagicMock()
        mock_mapper.expand_query.side_effect = RuntimeError("Expansion failed")
        mock_get_mapper.return_value = mock_mapper
        
        response = client.post(
            "/api/phenotypes/expand-query",
            json={"query": "Test query"},
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()

