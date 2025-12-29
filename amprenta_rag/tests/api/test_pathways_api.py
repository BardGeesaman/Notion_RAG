"""
Unit tests for pathways API endpoints.

Tests pathway enrichment and analysis with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestPathwayEnrich:
    """Tests for POST /api/pathways/enrich endpoint."""

    @patch("amprenta_rag.api.routers.pathways.perform_pathway_enrichment")
    @patch("amprenta_rag.api.routers.pathways._load_dataset_features")
    def test_pathway_enrich_success(self, mock_load_features, mock_enrichment):
        """Test successful pathway enrichment."""
        dataset_id = uuid4()
        
        # Mock features loading
        mock_load_features.return_value = (
            {"gene1", "gene2", "gene3"},
            {"gene"}
        )
        
        # Mock enrichment results
        mock_pathway = MagicMock()
        mock_pathway.pathway_id = "path:0001"
        mock_pathway.name = "Apoptosis Pathway"
        mock_pathway.source = "KEGG"
        
        mock_result = MagicMock()
        mock_result.pathway = mock_pathway
        mock_result.adjusted_p_value = 0.001
        mock_result.enrichment_ratio = 2.5
        mock_result.matched_features = ["gene1", "gene2"]
        
        mock_enrichment.return_value = [mock_result]
        
        response = client.post(
            f"/api/v1/pathways/enrich?dataset_id={dataset_id}",
        )
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["pathway_id"] == "path:0001"
        assert data[0]["name"] == "Apoptosis Pathway"
        assert data[0]["adjusted_p_value"] == 0.001

    @patch("amprenta_rag.api.routers.pathways._load_dataset_features")
    def test_pathway_enrich_no_features(self, mock_load_features):
        """Test pathway enrichment with no features."""
        dataset_id = uuid4()
        
        # Mock empty features
        mock_load_features.return_value = (set(), set())
        
        response = client.post(
            f"/api/v1/pathways/enrich?dataset_id={dataset_id}",
        )
        
        assert response.status_code == 404
        assert "no features" in response.json()["detail"].lower()


class TestProgramPathwayAnalysis:
    """Tests for GET /api/programs/{program_id}/pathway-analysis endpoint."""

    @patch("amprenta_rag.api.routers.pathways.get_cross_omics_enrichment")
    def test_program_pathway_analysis_success(self, mock_enrichment):
        """Test successful cross-omics pathway analysis."""
        program_id = uuid4()
        
        # Mock cross-omics result matching CrossOmicsEnrichmentResult schema
        mock_result = MagicMock()
        mock_result.asdict.return_value = {
            "program_id": str(program_id),
            "pathways": [
                {
                    "pathway_id": "path:0001",
                    "name": "Cell Cycle",
                    "source": "KEGG",
                    "adjusted_p_value": 0.0001,
                    "enrichment_ratio": 3.2,
                    "matched_by_omics": {
                        "transcriptomics": ["gene1", "gene2"],
                        "proteomics": ["protein1"],
                    },
                }
            ],
            "features": {
                "transcriptomics": ["gene1", "gene2"],
                "proteomics": ["protein1"],
                "metabolomics": [],
                "lipidomics": [],
            },
        }
        mock_enrichment.return_value = mock_result
        
        response = client.get(f"/api/v1/programs/{program_id}/pathway-analysis")
        
        assert response.status_code == 200
        data = response.json()
        assert data["program_id"] == str(program_id)
        assert len(data["pathways"]) == 1
        assert data["pathways"][0]["pathway_id"] == "path:0001"



