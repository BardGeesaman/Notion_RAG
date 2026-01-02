"""Tests for variant annotation API endpoints."""

import pytest
from uuid import uuid4
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient
from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


class TestAnnotationEndpoints:
    """Test annotation API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client with mocked auth."""
        mock_user = MagicMock()
        mock_user.id = uuid4()
        app.dependency_overrides[get_current_user] = lambda: mock_user
        try:
            yield TestClient(app)
        finally:
            app.dependency_overrides.clear()
    
    @patch("amprenta_rag.ingestion.genomics.annotation_service.annotate_variant_by_id")
    def test_annotate_single_variant(self, mock_annotate, client):
        """Test POST /variants/{variant_id}/annotate endpoint."""
        variant_id = uuid4()
        mock_annotation = MagicMock()
        mock_annotation.id = uuid4()
        mock_annotation.variant_id = variant_id
        mock_annotation.consequence = "missense_variant"
        mock_annotation.impact = "HIGH"
        mock_annotation.symbol = "TP53"
        mock_annotation.gene_id = "ENSG00000141510"
        mock_annotation.sift_prediction = "deleterious"
        mock_annotation.sift_score = 0.01
        mock_annotation.polyphen_prediction = "probably_damaging"
        mock_annotation.polyphen_score = 0.95
        mock_annotation.clin_sig = "pathogenic"
        mock_annotation.annotation_source = "VEP"
        mock_annotation.annotated_at = None
        mock_annotate.return_value = mock_annotation
        
        response = client.post(f"/api/v1/genomics/variants/{variant_id}/annotate")
        
        assert response.status_code == 200
        data = response.json()
        assert data["consequence"] == "missense_variant"
        assert data["impact"] == "HIGH"
    
    @patch("amprenta_rag.jobs.tasks.genomics.annotate_variants_task")
    def test_batch_annotation_job(self, mock_task, client):
        """Test POST /variants/annotate/batch endpoint."""
        mock_task.delay.return_value.id = "celery-job-123"
        
        variant_ids = [str(uuid4()) for _ in range(5)]
        
        response = client.post(
            "/api/v1/genomics/variants/annotate/batch",
            json={"variant_ids": variant_ids}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["job_id"] == "celery-job-123"
        assert data["queued_count"] == 5
    
    def test_batch_annotation_empty_list(self, client):
        """Test batch annotation with empty variant list."""
        response = client.post(
            "/api/v1/genomics/variants/annotate/batch",
            json={"variant_ids": []}
        )
        
        assert response.status_code == 400
    
    @patch("amprenta_rag.api.routers.genomics.db_session")
    def test_get_annotation_stats(self, mock_db_session, client):
        """Test GET /variants/annotations/stats endpoint."""
        mock_db = MagicMock()
        mock_db.query.return_value.scalar.return_value = 100
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/v1/genomics/variants/annotations/stats")
        
        assert response.status_code == 200
        data = response.json()
        assert "total_variants" in data