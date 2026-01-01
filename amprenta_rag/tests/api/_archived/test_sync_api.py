"""
Unit tests for sync API endpoints.

Tests external source sync operations with mocked dependencies.
"""

from __future__ import annotations

from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestSyncRun:
    """Tests for POST /api/sync/run endpoint."""

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_run_sync_success(self, mock_manager_class):
        """Test successful sync job creation."""
        job_id = uuid4()
        
        # Mock sync job
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.status = "pending"
        
        # Mock manager
        mock_manager = MagicMock()
        mock_manager.create_sync_job.return_value = mock_job
        mock_manager_class.return_value = mock_manager
        
        response = client.post(
            "/api/sync/run",
            json={"source": "chembl", "sync_type": "incremental"}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == str(job_id)
        assert data["status"] == "pending"

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_run_sync_invalid_source(self, mock_manager_class):
        """Test sync with invalid source."""
        # Mock manager to raise ValueError
        mock_manager = MagicMock()
        mock_manager.create_sync_job.side_effect = ValueError("Invalid source")
        mock_manager_class.return_value = mock_manager
        
        response = client.post(
            "/api/sync/run",
            json={"source": "invalid_source", "sync_type": "incremental"}
        )
        
        assert response.status_code == 400
        assert "invalid source" in response.json()["detail"].lower()


class TestSyncJobGet:
    """Tests for GET /api/sync/jobs/{job_id} endpoint."""

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_get_job_success(self, mock_manager_class):
        """Test successful job status retrieval."""
        job_id = uuid4()
        
        # Mock manager
        mock_manager = MagicMock()
        mock_manager.get_job_status.return_value = {
            "id": str(job_id),
            "source": "chembl",
            "status": "completed",
            "records_synced": 100,
            "records_new": 50,
            "records_updated": 50,
        }
        mock_manager_class.return_value = mock_manager
        
        response = client.get(f"/api/sync/jobs/{job_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(job_id)
        assert data["source"] == "chembl"
        assert data["status"] == "completed"
        assert data["records_synced"] == 100

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_get_job_not_found(self, mock_manager_class):
        """Test getting non-existent job."""
        job_id = uuid4()
        
        # Mock manager to raise ValueError
        mock_manager = MagicMock()
        mock_manager.get_job_status.side_effect = ValueError("Job not found")
        mock_manager_class.return_value = mock_manager
        
        response = client.get(f"/api/sync/jobs/{job_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestSyncJobsList:
    """Tests for GET /api/sync/jobs endpoint."""

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_list_jobs_success(self, mock_db_session):
        """Test successful jobs list retrieval."""
        job1_id = uuid4()
        job2_id = uuid4()
        
        # Mock jobs
        mock_job1 = MagicMock()
        mock_job1.id = job1_id
        mock_job1.source = "chembl"
        mock_job1.sync_type = "incremental"
        mock_job1.status = "completed"
        mock_job1.records_synced = 100
        mock_job1.records_updated = 50
        mock_job1.records_new = 50
        mock_job1.conflicts_detected = 0
        mock_job1.started_at = datetime(2024, 1, 1, 10, 0, 0, tzinfo=timezone.utc)
        mock_job1.completed_at = datetime(2024, 1, 1, 10, 30, 0, tzinfo=timezone.utc)
        mock_job1.created_at = datetime(2024, 1, 1, 9, 55, 0, tzinfo=timezone.utc)
        mock_job1.updated_at = datetime(2024, 1, 1, 10, 30, 0, tzinfo=timezone.utc)
        
        mock_job2 = MagicMock()
        mock_job2.id = job2_id
        mock_job2.source = "pubchem"
        mock_job2.sync_type = "full"
        mock_job2.status = "pending"
        mock_job2.records_synced = 0
        mock_job2.records_updated = 0
        mock_job2.records_new = 0
        mock_job2.conflicts_detected = 0
        mock_job2.started_at = None
        mock_job2.completed_at = None
        mock_job2.created_at = datetime(2024, 1, 2, 9, 0, 0, tzinfo=timezone.utc)
        mock_job2.updated_at = datetime(2024, 1, 2, 9, 0, 0, tzinfo=timezone.utc)
        
        # Mock database session
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.order_by.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.count.return_value = 2
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = [mock_job1, mock_job2]
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/sync/jobs?skip=0&limit=50")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 2
        assert len(data["jobs"]) == 2
        assert data["jobs"][0]["id"] == str(job1_id)
        assert data["jobs"][0]["source"] == "chembl"
        assert data["jobs"][1]["id"] == str(job2_id)
        assert data["jobs"][1]["source"] == "pubchem"

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_list_jobs_with_source_filter(self, mock_db_session):
        """Test jobs list with source filter."""
        job_id = uuid4()
        
        # Mock job
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.source = "chembl"
        mock_job.sync_type = "incremental"
        mock_job.status = "completed"
        mock_job.records_synced = 50
        mock_job.records_updated = 25
        mock_job.records_new = 25
        mock_job.conflicts_detected = 0
        mock_job.started_at = None
        mock_job.completed_at = None
        mock_job.created_at = None
        mock_job.updated_at = None
        
        # Mock database session
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.order_by.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.count.return_value = 1
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = [mock_job]
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/sync/jobs?source=chembl")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 1
        assert data["jobs"][0]["source"] == "chembl"


class TestSyncConflictsList:
    """Tests for GET /api/sync/conflicts endpoint."""

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_list_conflicts_success(self, mock_db_session):
        """Test successful conflicts list retrieval."""
        conflict_id = uuid4()
        record_id = uuid4()
        
        # Mock conflict and record
        mock_conflict = MagicMock()
        mock_conflict.id = conflict_id
        mock_conflict.record_id = record_id
        mock_conflict.conflict_type = "field_mismatch"
        mock_conflict.resolution_status = "pending"
        mock_conflict.resolved_at = None
        mock_conflict.local_value = "value_a"
        mock_conflict.external_value = "value_b"
        
        mock_record = MagicMock()
        mock_record.id = record_id
        mock_record.source = "chembl"
        mock_record.external_id = "CHEMBL12345"
        
        # Mock database session
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.join.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.count.return_value = 1
        mock_query.order_by.return_value = mock_query
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = [(mock_conflict, mock_record)]
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/sync/conflicts")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 1
        assert len(data["conflicts"]) == 1
        assert data["conflicts"][0]["id"] == str(conflict_id)
        assert data["conflicts"][0]["source"] == "chembl"
        assert data["conflicts"][0]["resolution_status"] == "pending"

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_list_conflicts_empty(self, mock_db_session):
        """Test listing conflicts when none exist."""
        # Mock database session with empty result
        mock_db = MagicMock()
        mock_query = mock_db.query.return_value
        mock_query.join.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.count.return_value = 0
        mock_query.order_by.return_value = mock_query
        mock_query.offset.return_value = mock_query
        mock_query.limit.return_value = mock_query
        mock_query.all.return_value = []
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.get("/api/sync/conflicts")
        
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 0
        assert data["conflicts"] == []


class TestSyncConflictResolve:
    """Tests for POST /api/sync/conflicts/{conflict_id}/resolve endpoint."""

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_resolve_conflict_success(self, mock_db_session):
        """Test successful conflict resolution."""
        conflict_id = uuid4()
        
        # Mock conflict
        mock_conflict = MagicMock()
        mock_conflict.id = conflict_id
        mock_conflict.resolution_status = "auto_merged"
        mock_conflict.resolved_at = datetime.now(timezone.utc)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_conflict
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.post(
            f"/api/sync/conflicts/{conflict_id}/resolve",
            json={"resolution": "auto_merged"}
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(conflict_id)
        assert data["resolution_status"] == "auto_merged"
        assert data["resolved_at"] is not None

    @patch("amprenta_rag.api.routers.sync.db_session")
    def test_resolve_conflict_not_found(self, mock_db_session):
        """Test resolving non-existent conflict."""
        conflict_id = uuid4()
        
        # Mock database session - conflict not found
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        response = client.post(
            f"/api/sync/conflicts/{conflict_id}/resolve",
            json={"resolution": "auto_merged"}
        )
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    def test_resolve_conflict_invalid_resolution(self):
        """Test resolving conflict with invalid resolution type."""
        conflict_id = uuid4()
        
        response = client.post(
            f"/api/sync/conflicts/{conflict_id}/resolve",
            json={"resolution": "invalid_type"}
        )
        
        # Should fail validation before hitting the database
        assert response.status_code == 422  # Validation error

