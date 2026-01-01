"""Integration tests for job queue management API with real database."""

import pytest
from datetime import datetime, timezone
from unittest.mock import patch
from uuid import uuid4

from amprenta_rag.database.models import PipelineJob, DockingRun, ExtractionJob, SyncJob, SingleCellDataset


@pytest.mark.integration
class TestJobsAPIIntegration:
    """Integration tests for job queue management endpoints."""

    def test_list_jobs_returns_empty(self, integration_client, db_session, timed_request):
        """Test that list_jobs returns empty list when no jobs exist."""
        response, benchmark = timed_request("GET", "/api/v1/jobs", "test_list_jobs_empty")
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["jobs"] == []
        assert data["total"] == 0

    def test_list_jobs_filters_by_type(self, integration_client, db_session, timed_request):
        """Test that job_type parameter filters results."""
        # Create real genomics job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc),
            progress_percent=100
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            "/api/v1/jobs?job_type=genomics", 
            "test_list_jobs_filter_type"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert len(data["jobs"]) >= 1
        found_job = next((j for j in data["jobs"] if j["id"] == str(job.id)), None)
        assert found_job is not None
        assert found_job["type"] == "genomics"
        assert found_job["status"] == "completed"

    def test_list_jobs_filters_by_status(self, integration_client, db_session, timed_request):
        """Test that status parameter filters results."""
        # Create jobs with different statuses
        running_job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="running",
            created_at=datetime.now(timezone.utc),
            progress_percent=50
        )
        completed_job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc),
            progress_percent=100
        )
        db_session.add_all([running_job, completed_job])
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET", 
            "/api/v1/jobs?status=running", 
            "test_list_jobs_filter_status"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        running_jobs = [j for j in data["jobs"] if j["status"] == "running"]
        assert len(running_jobs) >= 1
        found_job = next((j for j in running_jobs if j["id"] == str(running_job.id)), None)
        assert found_job is not None

    def test_get_job_success(self, integration_client, db_session, timed_request):
        """Test that get_job returns job details for valid ID."""
        # Create real job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc),
            progress_percent=100
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/genomics/{job.id}", 
            "test_get_job_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["id"] == str(job.id)
        assert data["type"] == "genomics"
        assert data["status"] == "completed"
        assert data["tool"] == "salmon"

    def test_get_job_not_found(self, integration_client, db_session, timed_request):
        """Test that get_job returns 404 for unknown job ID."""
        job_id = uuid4()
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/genomics/{job_id}", 
            "test_get_job_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "not found" in response.json()["detail"].lower()

    def test_cancel_job_success(self, integration_client, db_session, timed_request):
        """Test that cancel_job updates status to cancelled."""
        # Create running job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="running",
            created_at=datetime.now(timezone.utc),
            progress_percent=50
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/genomics/{job.id}/cancel", 
            "test_cancel_job_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "cancelled" in data["message"]
        assert data["status"] == "cancelled"
        
        # Verify job status was updated in database
        db_session.refresh(job)
        assert job.status == "cancelled"

    def test_cancel_pending_job(self, integration_client, db_session, timed_request):
        """Test that cancel_job works for pending status."""
        # Create pending job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="pending",
            created_at=datetime.now(timezone.utc),
            progress_percent=0
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/genomics/{job.id}/cancel", 
            "test_cancel_pending_job"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["status"] == "cancelled"
        
        # Verify real database state
        db_session.refresh(job)
        assert job.status == "cancelled"

    def test_cancel_completed_returns_error(self, integration_client, db_session, timed_request):
        """Test that cancel_job returns 400 for completed jobs."""
        # Create completed job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc),
            progress_percent=100
        )
        db_session.add(job)
        db_session.commit()
        
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/genomics/{job.id}/cancel", 
            "test_cancel_completed_error"
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "cannot cancel" in response.json()["detail"].lower()

    def test_cancel_not_found(self, integration_client, db_session, timed_request):
        """Test that cancel_job returns 404 for non-existent job."""
        job_id = uuid4()
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/docking/{job_id}/cancel", 
            "test_cancel_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "not found" in response.json()["detail"].lower()

    @patch('amprenta_rag.api.routers.jobs.celery_app')
    def test_retry_failed_job(self, mock_celery, integration_client, db_session, timed_request):
        """Test that retry_job attempts to resubmit failed job."""
        # Create failed job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="failed",
            created_at=datetime.now(timezone.utc),
            error_message="Some error"
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        # Mock Celery task submission (external service)
        mock_celery.send_task.return_value.id = "task-123"
        
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/genomics/{job.id}/retry", 
            "test_retry_failed_job"
        )
        
        # Should succeed with mocked Celery
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "resubmitted" in data["message"]
        
        # Verify job status was updated in database
        db_session.refresh(job)
        assert job.status in ["pending", "running", "queued"]

    def test_retry_non_failed_returns_error(self, integration_client, db_session, timed_request):
        """Test that retry_job returns 400 for non-failed jobs."""
        # Create completed job in database
        job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc),
            progress_percent=100
        )
        db_session.add(job)
        db_session.commit()
        
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/genomics/{job.id}/retry", 
            "test_retry_non_failed_error"
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "cannot retry" in response.json()["detail"].lower()

    def test_retry_not_found(self, integration_client, db_session, timed_request):
        """Test that retry_job returns 404 for non-existent job."""
        job_id = uuid4()
        response, benchmark = timed_request(
            "POST", 
            f"/api/v1/jobs/extraction/{job_id}/retry", 
            "test_retry_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "not found" in response.json()["detail"].lower()

    def test_invalid_job_type_returns_error(self, integration_client, timed_request):
        """Test that invalid job type returns 400 error."""
        job_id = uuid4()
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/invalid_type/{job_id}", 
            "test_invalid_job_type"
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "unknown job type" in response.json()["detail"].lower()

    def test_list_jobs_pagination(self, integration_client, db_session, timed_request):
        """Test that skip and limit parameters work correctly."""
        # Create multiple jobs
        jobs = []
        for i in range(15):
            job = PipelineJob(
                id=uuid4(),
                tool="salmon",
                status="completed",
                created_at=datetime.now(timezone.utc),
                progress_percent=100
            )
            jobs.append(job)
        
        db_session.add_all(jobs)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET", 
            "/api/v1/jobs?skip=10&limit=25", 
            "test_list_jobs_pagination"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["skip"] == 10
        assert data["limit"] == 25

    def test_get_job_docking_serialization(self, integration_client, db_session, timed_request):
        """Test docking job serialization includes docking-specific fields."""
        # Create real docking run in database
        job = DockingRun(
            id=uuid4(),
            status="completed",
            created_at=datetime.now(timezone.utc),
            started_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            total_compounds=100,
            completed_compounds=95,
            error_log="5 compounds failed"
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/docking/{job.id}", 
            "test_docking_serialization"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["type"] == "docking"
        assert data["total_compounds"] == 100
        assert data["completed_compounds"] == 95
        assert data["error_log"] == "5 compounds failed"

    def test_get_job_extraction_serialization(self, integration_client, db_session, timed_request):
        """Test extraction job serialization includes extraction-specific fields."""
        # Create real extraction job in database
        job = ExtractionJob(
            id=uuid4(),
            status="running",
            created_at=datetime.now(timezone.utc),
            started_at=datetime.now(timezone.utc),
            file_count=50,
            completed_count=30
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/extraction/{job.id}", 
            "test_extraction_serialization"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["type"] == "extraction"
        assert data["file_count"] == 50
        assert data["completed_count"] == 30

    def test_get_job_sync_serialization(self, integration_client, db_session, timed_request):
        """Test sync job serialization includes sync-specific fields."""
        # Create real sync job in database
        job = SyncJob(
            id=uuid4(),
            status="completed",
            created_at=datetime.now(timezone.utc),
            started_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            source="chembl",
            sync_type="incremental",
            records_synced=1000,
            records_new=50,
            records_updated=25
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/sync/{job.id}", 
            "test_sync_serialization"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["type"] == "sync"
        assert data["source"] == "chembl"
        assert data["sync_type"] == "incremental"
        assert data["records_synced"] == 1000
        assert data["records_new"] == 50
        assert data["records_updated"] == 25

    def test_get_job_single_cell_serialization(self, integration_client, db_session, timed_request):
        """Test single_cell job serialization includes single_cell-specific fields."""
        # Create real single cell dataset in database
        job = SingleCellDataset(
            id=uuid4(),
            processing_status="completed",
            created_at=datetime.now(timezone.utc),
            n_cells=5000,
            n_genes=20000,
            processing_log="QC complete, clustering done"
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        response, benchmark = timed_request(
            "GET", 
            f"/api/v1/jobs/single_cell/{job.id}", 
            "test_single_cell_serialization"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["type"] == "single_cell"
        assert data["n_cells"] == 5000
        assert data["n_genes"] == 20000
        assert data["processing_log"] == "QC complete, clustering done"

    def test_list_jobs_multiple_types_combined(self, integration_client, db_session, timed_request):
        """Test listing all job types without filter returns combined results."""
        # Create jobs of different types
        genomics_job = PipelineJob(
            id=uuid4(),
            tool="salmon",
            status="completed",
            created_at=datetime.now(timezone.utc)
        )
        extraction_job = ExtractionJob(
            id=uuid4(),
            status="running",
            created_at=datetime.now(timezone.utc),
            file_count=10
        )
        sync_job = SyncJob(
            id=uuid4(),
            status="pending",
            created_at=datetime.now(timezone.utc),
            source="chembl"
        )
        
        db_session.add_all([genomics_job, extraction_job, sync_job])
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET", 
            "/api/v1/jobs", 
            "test_list_multiple_types"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        job_types = {job["type"] for job in data["jobs"]}
        # Should have at least the types we created
        assert "genomics" in job_types
        assert "extraction" in job_types
        assert "sync" in job_types
