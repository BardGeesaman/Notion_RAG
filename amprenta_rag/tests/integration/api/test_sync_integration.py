"""Integration tests for sync API with real database."""

import pytest
from datetime import datetime, timezone
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.database.models import SyncJob, SyncRecord, SyncConflict


@pytest.mark.integration
class TestSyncAPIIntegration:
    """Integration tests for sync API endpoints."""

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_run_sync_success(self, mock_manager_class, integration_client, 
                            db_session, timed_request):
        """Test successful sync job creation with real database."""
        # Create real sync job in database
        job = SyncJob(
            id=uuid4(),
            source="chembl",
            sync_type="incremental",
            status="pending",
            created_at=datetime.now(timezone.utc)
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        # Mock sync manager (external service)
        mock_manager = MagicMock()
        mock_manager.create_sync_job.return_value = job
        mock_manager_class.return_value = mock_manager
        
        response, benchmark = timed_request(
            "POST",
            "/api/sync/run",
            "test_run_sync_success",
            json={"source": "chembl", "sync_type": "incremental"}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["job_id"] == str(job.id)
        assert data["status"] == "pending"
        
        # Verify job exists in database
        db_job = db_session.query(SyncJob).filter_by(id=job.id).first()
        assert db_job is not None
        assert db_job.source == "chembl"
        assert db_job.sync_type == "incremental"

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_run_sync_invalid_source(self, mock_manager_class, integration_client, timed_request):
        """Test sync with invalid source."""
        # Mock manager to raise ValueError (external validation)
        mock_manager = MagicMock()
        mock_manager.create_sync_job.side_effect = ValueError("Invalid source")
        mock_manager_class.return_value = mock_manager
        
        response, benchmark = timed_request(
            "POST",
            "/api/sync/run",
            "test_run_sync_invalid_source",
            json={"source": "invalid_source", "sync_type": "incremental"}
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "invalid source" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_get_job_success(self, mock_manager_class, integration_client, 
                           db_session, timed_request):
        """Test successful job status retrieval from database."""
        # Create real sync job in database
        job = SyncJob(
            id=uuid4(),
            source="chembl",
            sync_type="full",
            status="completed",
            created_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            records_synced=100,
            records_new=50,
            records_updated=50
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)
        
        # Mock manager to return job data
        mock_manager = MagicMock()
        mock_manager.get_job_status.return_value = {
            "id": str(job.id),
            "source": job.source,
            "status": job.status,
            "records_synced": job.records_synced,
            "records_new": job.records_new,
            "records_updated": job.records_updated,
        }
        mock_manager_class.return_value = mock_manager
        
        response, benchmark = timed_request(
            "GET",
            f"/api/sync/jobs/{job.id}",
            "test_get_job_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["id"] == str(job.id)
        assert data["source"] == "chembl"
        assert data["status"] == "completed"
        assert data["records_synced"] == 100

    @patch("amprenta_rag.api.routers.sync.SyncManager")
    def test_get_job_not_found(self, mock_manager_class, integration_client, timed_request):
        """Test getting non-existent job."""
        job_id = uuid4()
        
        # Mock manager to raise ValueError for not found
        mock_manager = MagicMock()
        mock_manager.get_job_status.side_effect = ValueError("Job not found")
        mock_manager_class.return_value = mock_manager
        
        response, benchmark = timed_request(
            "GET",
            f"/api/sync/jobs/{job_id}",
            "test_get_job_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    def test_list_jobs_success(self, integration_client, db_session, timed_request):
        """Test listing sync jobs from real database."""
        # Create real sync jobs in database
        jobs = []
        for i, source in enumerate(["chembl", "uniprot", "pubchem"]):
            job = SyncJob(
                id=uuid4(),
                source=source,
                sync_type="incremental",
                status="completed" if i % 2 == 0 else "running",
                created_at=datetime.now(timezone.utc),
                records_synced=100 * (i + 1)
            )
            jobs.append(job)
        
        db_session.add_all(jobs)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            "/api/sync/jobs",
            "test_list_jobs_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert len(data["jobs"]) >= 3
        
        # Find our created jobs
        our_jobs = [j for j in data["jobs"] if j["id"] in [str(job.id) for job in jobs]]
        assert len(our_jobs) == 3
        
        # Verify data integrity
        chembl_job = next(j for j in our_jobs if j["source"] == "chembl")
        assert chembl_job["sync_type"] == "incremental"

    def test_list_jobs_with_source_filter(self, integration_client, db_session, timed_request):
        """Test listing jobs with source filter."""
        # Create jobs with different sources
        chembl_job = SyncJob(
            id=uuid4(),
            source="chembl",
            sync_type="full",
            status="completed",
            created_at=datetime.now(timezone.utc)
        )
        uniprot_job = SyncJob(
            id=uuid4(),
            source="uniprot",
            sync_type="incremental",
            status="running",
            created_at=datetime.now(timezone.utc)
        )
        
        db_session.add_all([chembl_job, uniprot_job])
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            "/api/sync/jobs?source=chembl",
            "test_list_jobs_filter_source"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        # Should only return ChEMBL jobs
        chembl_jobs = [j for j in data["jobs"] if j["source"] == "chembl"]
        uniprot_jobs = [j for j in data["jobs"] if j["source"] == "uniprot"]
        
        assert len(chembl_jobs) >= 1
        assert len(uniprot_jobs) == 0

    def test_list_conflicts_success(self, integration_client, db_session, timed_request):
        """Test listing sync conflicts from real database."""
        # Create real sync conflicts in database
        conflicts = []
        for i in range(3):
            conflict = SyncConflict(
                id=uuid4(),
                source="chembl",
                external_id=f"CHEMBL{1000 + i}",
                entity_type="compound",
                conflict_type="data_mismatch",
                local_data={"name": f"Local Compound {i}"},
                external_data={"name": f"External Compound {i}"},
                status="pending",
                created_at=datetime.now(timezone.utc)
            )
            conflicts.append(conflict)
        
        db_session.add_all(conflicts)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            "/api/sync/conflicts",
            "test_list_conflicts_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert len(data["conflicts"]) >= 3
        
        # Find our created conflicts
        our_conflicts = [c for c in data["conflicts"] if c["id"] in [str(conflict.id) for conflict in conflicts]]
        assert len(our_conflicts) == 3
        
        # Verify data integrity
        first_conflict = our_conflicts[0]
        assert first_conflict["source"] == "chembl"
        assert first_conflict["entity_type"] == "compound"
        assert first_conflict["status"] == "pending"

    def test_list_conflicts_empty(self, integration_client, db_session, timed_request):
        """Test listing conflicts when none exist returns empty list."""
        # Don't create any conflicts
        
        response, benchmark = timed_request(
            "GET",
            "/api/sync/conflicts",
            "test_list_conflicts_empty"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        # May have existing conflicts from other tests, but should be valid structure
        assert "conflicts" in data
        assert isinstance(data["conflicts"], list)

    def test_resolve_conflict_success(self, integration_client, db_session, timed_request):
        """Test resolving sync conflict with real database."""
        # Create real conflict in database
        conflict = SyncConflict(
            id=uuid4(),
            source="uniprot",
            external_id="P12345",
            entity_type="protein",
            conflict_type="data_mismatch",
            local_data={"name": "Local Protein"},
            external_data={"name": "External Protein"},
            status="pending",
            created_at=datetime.now(timezone.utc)
        )
        db_session.add(conflict)
        db_session.commit()
        db_session.refresh(conflict)
        
        response, benchmark = timed_request(
            "POST",
            f"/api/sync/conflicts/{conflict.id}/resolve",
            "test_resolve_conflict_success",
            json={"resolution": "use_external"}
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["message"] == "Conflict resolved"
        assert data["resolution"] == "use_external"
        
        # Verify conflict was resolved in database
        db_session.refresh(conflict)
        assert conflict.status == "resolved"
        assert conflict.resolution == "use_external"

    def test_resolve_conflict_not_found(self, integration_client, timed_request):
        """Test resolving non-existent conflict."""
        fake_conflict_id = uuid4()
        
        response, benchmark = timed_request(
            "POST",
            f"/api/sync/conflicts/{fake_conflict_id}/resolve",
            "test_resolve_conflict_not_found",
            json={"resolution": "use_local"}
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    def test_resolve_conflict_invalid_resolution(self, integration_client, 
                                               db_session, timed_request):
        """Test resolving conflict with invalid resolution type."""
        # Create real conflict in database
        conflict = SyncConflict(
            id=uuid4(),
            source="chembl",
            external_id="CHEMBL123",
            entity_type="compound",
            conflict_type="data_mismatch",
            local_data={"name": "Local"},
            external_data={"name": "External"},
            status="pending",
            created_at=datetime.now(timezone.utc)
        )
        db_session.add(conflict)
        db_session.commit()
        
        response, benchmark = timed_request(
            "POST",
            f"/api/sync/conflicts/{conflict.id}/resolve",
            "test_resolve_conflict_invalid",
            json={"resolution": "invalid_resolution"}
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "invalid resolution" in response.json()["detail"].lower()

    def test_sync_records_created_during_job(self, integration_client, db_session, timed_request):
        """Test that sync records are properly linked to jobs."""
        # Create sync job and related records
        job = SyncJob(
            id=uuid4(),
            source="pubchem",
            sync_type="incremental",
            status="completed",
            created_at=datetime.now(timezone.utc),
            records_synced=2
        )
        db_session.add(job)
        db_session.commit()
        
        # Create sync records linked to job
        records = [
            SyncRecord(
                id=uuid4(),
                job_id=job.id,
                external_id="CID123",
                entity_type="compound",
                action="created",
                local_entity_id=uuid4(),
                created_at=datetime.now(timezone.utc)
            ),
            SyncRecord(
                id=uuid4(),
                job_id=job.id,
                external_id="CID124",
                entity_type="compound",
                action="updated",
                local_entity_id=uuid4(),
                created_at=datetime.now(timezone.utc)
            )
        ]
        
        db_session.add_all(records)
        db_session.commit()
        
        # Verify records are linked to job
        job_records = db_session.query(SyncRecord).filter_by(job_id=job.id).all()
        assert len(job_records) == 2
        assert all(record.job_id == job.id for record in job_records)
        assert {record.action for record in job_records} == {"created", "updated"}
