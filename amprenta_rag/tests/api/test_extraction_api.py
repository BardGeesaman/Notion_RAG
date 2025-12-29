"""
Unit tests for extraction API endpoints.

Tests batch extraction job operations with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestExtractionJobGet:
    """Tests for GET /api/extraction/jobs/{job_id} endpoint."""

    @patch("amprenta_rag.api.routers.extraction.ExtractionBatchService")
    def test_get_job_success(self, mock_service_class):
        """Test successful job status retrieval."""
        job_id = uuid4()
        
        # Mock service
        mock_service = MagicMock()
        mock_service.get_job_status.return_value = {
            "id": str(job_id),
            "status": "processing",
            "file_count": 5,
            "completed_count": 3,
            "progress_pct": 60.0,
        }
        mock_service_class.return_value = mock_service
        
        response = client.get(f"/api/extraction/jobs/{job_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(job_id)
        assert data["status"] == "processing"
        assert data["progress_pct"] == 60.0

    @patch("amprenta_rag.api.routers.extraction.ExtractionBatchService")
    def test_get_job_not_found(self, mock_service_class):
        """Test getting non-existent job."""
        job_id = uuid4()
        
        # Mock service to raise ValueError (job not found)
        mock_service = MagicMock()
        mock_service.get_job_status.side_effect = ValueError("Job not found")
        mock_service_class.return_value = mock_service
        
        response = client.get(f"/api/extraction/jobs/{job_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestExtractionJobsList:
    """Tests for GET /api/extraction/jobs endpoint."""

    @patch("amprenta_rag.api.routers.extraction.db_session")
    def test_list_jobs_success(self, mock_db_session):
        """Test successful jobs list retrieval."""
        job1_id = uuid4()
        job2_id = uuid4()
        
        # Mock jobs
        mock_job1 = MagicMock()
        mock_job1.id = job1_id
        mock_job1.batch_id = None
        mock_job1.file_count = 5
        mock_job1.completed_count = 5
        mock_job1.status = "completed"
        mock_job1.created_at = None
        mock_job1.updated_at = None
        
        mock_job2 = MagicMock()
        mock_job2.id = job2_id
        mock_job2.batch_id = None
        mock_job2.file_count = 10
        mock_job2.completed_count = 7
        mock_job2.status = "processing"
        mock_job2.created_at = None
        mock_job2.updated_at = None
        
        # Mock database session with multiple query chains
        mock_db = MagicMock()
        
        # First query chain: get jobs
        mock_jobs_query = MagicMock()
        mock_jobs_query.order_by.return_value = mock_jobs_query
        mock_jobs_query.count.return_value = 2
        mock_jobs_query.offset.return_value = mock_jobs_query
        mock_jobs_query.limit.return_value = mock_jobs_query
        mock_jobs_query.all.return_value = [mock_job1, mock_job2]
        
        # Second query chain: document counts (returns tuples of (job_id, count))
        mock_docs_query = MagicMock()
        mock_docs_query.filter.return_value = mock_docs_query
        mock_docs_query.group_by.return_value = mock_docs_query
        mock_docs_query.all.return_value = [(job1_id, 10), (job2_id, 15)]
        
        # db.query() returns different mocks based on what's queried
        def mock_query_side_effect(*args, **kwargs):
            # First arg is the model or models
            if args and "ExtractionJob" in str(args[0]):
                return mock_jobs_query
            else:
                return mock_docs_query
        
        mock_db.query.side_effect = mock_query_side_effect
        
        # Mock func.count for document count aggregation
        mock_db.func.count.return_value = MagicMock()
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/extraction/jobs?skip=0&limit=50")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 2
        assert len(data["jobs"]) == 2
        assert data["jobs"][0]["id"] == str(job1_id)
        assert data["jobs"][0]["status"] == "completed"
        assert data["jobs"][1]["id"] == str(job2_id)
        assert data["jobs"][1]["status"] == "processing"

    @patch("amprenta_rag.api.routers.extraction.db_session")
    def test_list_jobs_empty(self, mock_db_session):
        """Test listing jobs when none exist."""
        # Mock database session with empty result
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.order_by.return_value = mock_query
        mock_query.count.return_value = 0
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = []
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/extraction/jobs")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 0
        assert data["jobs"] == []

    @patch("amprenta_rag.api.routers.extraction.db_session")
    def test_list_jobs_with_pagination(self, mock_db_session):
        """Test jobs list with pagination parameters."""
        # Mock database session
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.order_by.return_value = mock_query
        mock_query.count.return_value = 100
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = []
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/extraction/jobs?skip=20&limit=10")
        
        assert response.status_code == 200
        data = response.json()
        assert data["skip"] == 20
        assert data["limit"] == 10
        assert data["total"] == 100

