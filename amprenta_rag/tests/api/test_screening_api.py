"""
Unit tests for screening API endpoints.

Tests HTS screening campaigns and active learning with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestListCampaigns:
    """Tests for GET /api/screening/campaigns endpoint."""

    @patch("amprenta_rag.api.services.screening.list_campaigns")
    def test_list_campaigns_success(self, mock_list):
        """Test successful campaigns list retrieval."""
        mock_list.return_value = [
            {
                "id": str(uuid4()),
                "campaign_id": "CAMP001",
                "campaign_name": "ALS Target Screen",
                "assay_type": "cell_viability",
                "target": "SOD1",
                "total_wells": 384,
                "hit_count": 15,
                "run_date": "2024-01-15",
            },
            {
                "id": str(uuid4()),
                "campaign_id": "CAMP002",
                "campaign_name": "Kinase Inhibitor Screen",
                "assay_type": "enzyme_activity",
                "target": "CDK2",
                "total_wells": 1536,
                "hit_count": 42,
                "run_date": "2024-02-01",
            },
        ]
        
        response = client.get("/api/v1/screening/campaigns")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["campaign_id"] == "CAMP001"
        assert data[1]["campaign_id"] == "CAMP002"

    @patch("amprenta_rag.api.services.screening.list_campaigns")
    def test_list_campaigns_empty(self, mock_list):
        """Test listing campaigns when none exist."""
        mock_list.return_value = []
        
        response = client.get("/api/v1/screening/campaigns")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []


class TestGetCampaign:
    """Tests for GET /api/screening/campaigns/{campaign_id} endpoint."""

    @patch("amprenta_rag.api.services.screening.get_campaign")
    def test_get_campaign_success(self, mock_get):
        """Test successful campaign retrieval."""
        campaign_id = "CAMP001"
        
        mock_get.return_value = {
            "id": str(uuid4()),
            "campaign_id": campaign_id,
            "campaign_name": "Test Campaign",
            "assay_type": "binding_assay",
            "target": "EGFR",
            "total_wells": 96,
            "hit_count": 5,
            "run_date": "2024-01-20",
        }
        
        response = client.get(f"/api/v1/screening/campaigns/{campaign_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["campaign_id"] == campaign_id
        assert data["campaign_name"] == "Test Campaign"

    @patch("amprenta_rag.api.services.screening.get_campaign")
    def test_get_campaign_not_found(self, mock_get):
        """Test getting non-existent campaign."""
        campaign_id = "NONEXISTENT"
        
        mock_get.return_value = None
        
        response = client.get(f"/api/v1/screening/campaigns/{campaign_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestGetCampaignHits:
    """Tests for GET /api/screening/campaigns/{campaign_id}/hits endpoint."""

    @patch("amprenta_rag.api.services.screening.get_campaign_hits")
    def test_get_campaign_hits_success(self, mock_get_hits):
        """Test successful campaign hits retrieval."""
        campaign_id = "CAMP001"
        
        mock_get_hits.return_value = [
            {
                "id": str(uuid4()),
                "result_id": "RES001",
                "compound_id": str(uuid4()),
                "well_position": "A1",
                "raw_value": 0.85,
                "normalized_value": 0.92,
                "z_score": 3.5,
                "hit_flag": True,
                "hit_category": "strong",
            },
            {
                "id": str(uuid4()),
                "result_id": "RES002",
                "compound_id": str(uuid4()),
                "well_position": "B3",
                "raw_value": 0.75,
                "normalized_value": 0.82,
                "z_score": 2.8,
                "hit_flag": True,
                "hit_category": "moderate",
            },
        ]
        
        response = client.get(f"/api/v1/screening/campaigns/{campaign_id}/hits")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["hit_category"] == "strong"
        assert data[1]["hit_category"] == "moderate"

    @patch("amprenta_rag.api.services.screening.get_campaign_hits")
    def test_get_campaign_hits_empty(self, mock_get_hits):
        """Test getting hits for campaign with no hits."""
        campaign_id = "CAMP001"
        
        mock_get_hits.return_value = []
        
        response = client.get(f"/api/v1/screening/campaigns/{campaign_id}/hits")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []


class TestSuggestCompounds:
    """Tests for POST /api/screening/suggest endpoint."""

    @patch("amprenta_rag.analysis.active_learning.suggest_next_compounds")
    def test_suggest_compounds_success(self, mock_suggest):
        """Test successful compound suggestions."""
        cmpd1_id = uuid4()
        cmpd2_id = uuid4()
        screened_id = uuid4()
        
        # Mock suggestions
        mock_suggestion1 = MagicMock()
        mock_suggestion1.compound_id = cmpd1_id
        mock_suggestion1.smiles = "CCO"
        mock_suggestion1.acquisition_score = 0.95
        mock_suggestion1.strategy_used = "uncertainty"
        mock_suggestion1.rank = 1
        mock_suggestion1.explanation = "High uncertainty"
        
        mock_suggestion2 = MagicMock()
        mock_suggestion2.compound_id = cmpd2_id
        mock_suggestion2.smiles = "CC(C)O"
        mock_suggestion2.acquisition_score = 0.88
        mock_suggestion2.strategy_used = "uncertainty"
        mock_suggestion2.rank = 2
        mock_suggestion2.explanation = "Moderate uncertainty"
        
        mock_suggest.return_value = [mock_suggestion1, mock_suggestion2]
        
        response = client.post(
            "/api/v1/screening/suggest",
            json={
                "screened": [
                    {"compound_id": str(screened_id), "smiles": "C", "activity": True}
                ],
                "candidates": [
                    {"compound_id": str(cmpd1_id), "smiles": "CCO"},
                    {"compound_id": str(cmpd2_id), "smiles": "CC(C)O"},
                ],
                "strategy": "uncertainty",
                "batch_size": 10,
                "model_id": None,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert len(data["suggestions"]) == 2
        assert data["suggestions"][0]["compound_id"] == str(cmpd1_id)
        assert data["suggestions"][1]["compound_id"] == str(cmpd2_id)
        assert data["strategy_used"] == "uncertainty"

    @patch("amprenta_rag.analysis.active_learning.suggest_next_compounds")
    def test_suggest_compounds_error(self, mock_suggest):
        """Test compound suggestions with error."""
        cmpd_id = uuid4()
        screened_id = uuid4()
        
        mock_suggest.side_effect = RuntimeError("Suggestion failed")
        
        response = client.post(
            "/api/v1/screening/suggest",
            json={
                "screened": [{"compound_id": str(screened_id), "smiles": "C", "activity": True}],
                "candidates": [{"compound_id": str(cmpd_id), "smiles": "CCO"}],
                "strategy": "diversity",
                "batch_size": 5,
                "model_id": None,
            },
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()


class TestSuggestCompoundsForProgram:
    """Tests for POST /api/screening/suggest/program/{program_id} endpoint."""

    @patch("amprenta_rag.analysis.active_learning.get_compound_suggestions_for_program")
    def test_suggest_for_program_success(self, mock_suggest):
        """Test successful program-based suggestions."""
        program_id = uuid4()
        cmpd_id = uuid4()
        
        # Mock suggestions
        mock_suggestion = MagicMock()
        mock_suggestion.compound_id = cmpd_id
        mock_suggestion.smiles = "CCO"
        mock_suggestion.acquisition_score = 0.92
        mock_suggestion.strategy_used = "diversity"
        mock_suggestion.rank = 1
        mock_suggestion.explanation = "Diverse structure"
        
        mock_suggest.return_value = [mock_suggestion]
        
        response = client.post(
            f"/api/v1/screening/suggest/program/{program_id}",
            json={
                "program_id": str(program_id),
                "strategy": "diversity",
                "batch_size": 5,
                "model_id": None,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert len(data["suggestions"]) == 1
        assert data["suggestions"][0]["compound_id"] == str(cmpd_id)
        assert data["strategy_used"] == "diversity"

    @patch("amprenta_rag.analysis.active_learning.get_compound_suggestions_for_program")
    def test_suggest_for_program_error(self, mock_suggest):
        """Test program-based suggestions with error."""
        program_id = uuid4()
        
        mock_suggest.side_effect = RuntimeError("No screening data for program")
        
        response = client.post(
            f"/api/v1/screening/suggest/program/{program_id}",
            json={
                "program_id": str(program_id),
                "strategy": "uncertainty",
                "batch_size": 10,
                "model_id": None,
            },
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()

