"""
Unit tests for batch correction API endpoints.

Tests batch effect correction with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestBatchCorrect:
    """Tests for POST /api/analysis/batch-correct endpoint."""

    @patch("amprenta_rag.api.routers.batch.correct_batch_effects")
    def test_batch_correct_success(self, mock_correct):
        """Test successful batch correction."""
        dataset1_id = uuid4()
        dataset2_id = uuid4()
        corrected_dataset_id = uuid4()
        
        # Mock correction result
        import pandas as pd
        mock_df = pd.DataFrame({
            "feature": ["gene1", "gene2", "gene3"],
            "sample1": [1.0, 2.0, 3.0],
            "sample2": [1.5, 2.5, 3.5],
        })
        
        mock_correct.return_value = {
            "corrected_df": mock_df,
            "corrected_dataset_id": corrected_dataset_id,
            "stats": {
                "method": "combat",
                "num_batches": 2,
                "num_features": 3,
                "variance_explained": 0.85,
            },
        }
        
        response = client.post(
            "/api/analysis/batch-correct",
            json={
                "datasets": [str(dataset1_id), str(dataset2_id)],
                "batch_map": {"batch1": [str(dataset1_id)], "batch2": [str(dataset2_id)]},
                "method": "combat",
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["corrected_dataset_id"] == str(corrected_dataset_id)
        assert data["stats"]["method"] == "combat"
        assert len(data["preview"]) > 0

    @patch("amprenta_rag.api.routers.batch.correct_batch_effects")
    def test_batch_correct_import_error(self, mock_correct):
        """Test batch correction with missing dependency."""
        mock_correct.side_effect = ImportError("combat module not available")
        
        response = client.post(
            "/api/analysis/batch-correct",
            json={
                "datasets": [str(uuid4())],
                "batch_map": {"batch1": [str(uuid4())]},
                "method": "combat",
            },
        )
        
        assert response.status_code == 503
        assert "combat" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.batch.correct_batch_effects")
    def test_batch_correct_value_error(self, mock_correct):
        """Test batch correction with invalid input."""
        mock_correct.side_effect = ValueError("Invalid batch configuration")
        
        response = client.post(
            "/api/analysis/batch-correct",
            json={
                "datasets": [str(uuid4())],
                "batch_map": {},
                "method": "combat",
            },
        )
        
        assert response.status_code == 422
        assert "invalid" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.batch.correct_batch_effects")
    def test_batch_correct_generic_error(self, mock_correct):
        """Test batch correction with generic error."""
        mock_correct.side_effect = RuntimeError("Unexpected error")
        
        response = client.post(
            "/api/analysis/batch-correct",
            json={
                "datasets": [str(uuid4())],
                "batch_map": {"batch1": [str(uuid4())]},
                "method": "combat",
            },
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()

