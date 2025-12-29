"""
Unit tests for predictors API endpoints.

Tests assay outcome prediction model training and inference with mocked dependencies.
"""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestTrainAssayPredictor:
    """Tests for POST /api/v1/predictors/train endpoint."""

    @patch("amprenta_rag.analysis.assay_predictor.train_assay_predictor")
    def test_train_predictor_success(self, mock_train):
        """Test successful predictor training."""
        program_id = uuid4()
        model_id = uuid4()
        
        # Mock training result
        mock_result = MagicMock()
        mock_result.model_id = model_id
        mock_result.program_id = program_id
        mock_result.assay_type = "cell_viability"
        mock_result.model_performance = {
            "accuracy": 0.85,
            "roc_auc": 0.88,
            "f1_score": 0.82,
        }
        mock_result.training_stats = MagicMock()
        mock_result.training_stats.total_compounds = 500
        mock_result.training_stats.active_compounds = 50
        mock_result.training_stats.inactive_compounds = 450
        mock_result.training_stats.activity_rate = 0.1
        mock_result.training_stats.feature_count = 2048
        mock_result.training_stats.data_quality_score = 0.90
        mock_result.feature_names = ["morgan_fp"]
        mock_result.training_time_seconds = 45.2
        mock_result.success = True
        mock_result.error_message = None
        
        mock_train.return_value = mock_result
        
        response = client.post(
            "/api/v1/predictors/train",
            json={
                "program_id": str(program_id),
                "assay_type": "cell_viability",
                "features": ["Morgan"],
                "min_actives": 50,
                "min_inactives": 50,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["model_id"] == str(model_id)
        assert data["success"] is True
        assert data["model_performance"]["accuracy"] == 0.85

    @patch("amprenta_rag.analysis.assay_predictor.train_assay_predictor")
    def test_train_predictor_error(self, mock_train):
        """Test predictor training with error."""
        program_id = uuid4()
        
        mock_train.side_effect = RuntimeError("Insufficient training data")
        
        response = client.post(
            "/api/v1/predictors/train",
            json={
                "program_id": str(program_id),
                "assay_type": "binding_assay",
                "features": ["Morgan"],
                "min_actives": 50,
                "min_inactives": 50,
            },
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()


class TestPredictAssayOutcome:
    """Tests for POST /api/v1/predictors/{model_id}/predict endpoint."""

    def test_predict_with_valid_smiles(self):
        """Test prediction endpoint accepts valid SMILES."""
        model_id = uuid4()
        
        # Endpoint returns 200 even if model doesn't exist (logs error, returns empty)
        response = client.post(
            f"/api/v1/predictors/{model_id}/predict",
            json={"compound_smiles": ["CCO", "CC(C)O"]},
        )
        
        # Real implementation returns 200 with empty predictions if model not found
        assert response.status_code == 200
        data = response.json()
        assert "predictions" in data
        assert "model_id" in data

    def test_predict_invalid_request(self):
        """Test prediction with invalid request."""
        model_id = uuid4()
        
        # Empty SMILES list should fail validation
        response = client.post(
            f"/api/v1/predictors/{model_id}/predict",
            json={"compound_smiles": []},
        )
        
        # Should get validation error or 500
        assert response.status_code in [422, 500]


class TestListAssayModels:
    """Tests for GET /api/v1/predictors endpoint."""

    def test_list_models_success(self):
        """Test successful model listing."""
        # Test with real implementation (no models will exist, returns empty)
        response = client.get("/api/v1/predictors")
        
        assert response.status_code == 200
        data = response.json()
        assert "total_models" in data
        assert "models" in data
        assert isinstance(data["models"], list)

    def test_list_models_with_program_filter(self):
        """Test model listing with program filter."""
        program_id = uuid4()
        
        # Test with real implementation
        response = client.get(f"/api/v1/predictors?program_id={program_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert "total_models" in data
        assert "models" in data

