"""Integration tests for backup API with real database."""

import pytest
import tempfile
from datetime import datetime, timezone, timedelta
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.database.models import BackupRecord, ProjectExport


@pytest.mark.integration
class TestBackupAPIIntegration:
    """Integration tests for backup API endpoints."""

    @patch('amprenta_rag.api.routers.backup.run_database_backup')
    def test_trigger_manual_backup_success(self, mock_backup_task, integration_client, 
                                         db_session, timed_request):
        """Test successful manual backup trigger with real database."""
        # Mock Celery task (external service)
        mock_task = MagicMock()
        mock_task.id = "task_12345"
        mock_backup_task.delay.return_value = mock_task
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/database",
            "test_trigger_backup",
            json={"backup_type": "full"}
        )
        
        assert response.status_code == 202
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["message"] == "Full backup initiated"
        assert data["task_id"] == "task_12345"
        assert data["backup_type"] == "full"
        
        # Verify task was called
        mock_backup_task.delay.assert_called_once_with("full")

    def test_trigger_manual_backup_invalid_type(self, integration_client, timed_request):
        """Test manual backup with invalid backup type."""
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/database",
            "test_trigger_backup_invalid",
            json={"backup_type": "invalid"}
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "Invalid backup type" in response.json()["detail"]

    def test_get_backup_history_success(self, integration_client, db_session, timed_request):
        """Test retrieving backup history from real database."""
        # Create real backup records in database
        backups = []
        for i, status in enumerate(["completed", "running", "failed"]):
            backup = BackupRecord(
                id=uuid4(),
                backup_type="full",
                status=status,
                created_at=datetime.now(timezone.utc) - timedelta(days=i),
                file_size=1000000 + i * 100000 if status == "completed" else None,
                file_path=f"backups/backup_{i}.sql" if status == "completed" else None
            )
            backups.append(backup)
        
        db_session.add_all(backups)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            "/api/v1/backup/history",
            "test_get_backup_history"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert len(data["backups"]) >= 3
        
        # Find our created backups
        our_backups = [b for b in data["backups"] if b["id"] in [str(backup.id) for backup in backups]]
        assert len(our_backups) == 3
        
        # Verify data integrity
        completed_backup = next(b for b in our_backups if b["status"] == "completed")
        assert completed_backup["backup_type"] == "full"
        assert completed_backup["file_size"] is not None

    def test_get_backup_details_success(self, integration_client, db_session, timed_request):
        """Test retrieving specific backup details from database."""
        # Create real backup record in database
        backup = BackupRecord(
            id=uuid4(),
            backup_type="incremental",
            status="completed",
            created_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            file_size=2500000,
            file_path="backups/backup_detailed.sql",
            metadata={"tables": 15, "rows": 50000}
        )
        db_session.add(backup)
        db_session.commit()
        db_session.refresh(backup)
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/{backup.id}",
            "test_get_backup_details"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["id"] == str(backup.id)
        assert data["backup_type"] == "incremental"
        assert data["status"] == "completed"
        assert data["file_size"] == 2500000
        assert data["metadata"]["tables"] == 15

    def test_get_backup_details_not_found(self, integration_client, timed_request):
        """Test retrieving non-existent backup returns 404."""
        fake_backup_id = uuid4()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/{fake_backup_id}",
            "test_get_backup_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    @patch('amprenta_rag.api.routers.backup.S3Client')
    def test_download_backup_success(self, mock_s3_client, integration_client, 
                                   db_session, timed_request):
        """Test downloading completed backup with mocked S3."""
        # Create completed backup in database
        backup = BackupRecord(
            id=uuid4(),
            backup_type="full",
            status="completed",
            created_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            file_size=1000000,
            file_path="backups/download_test.sql"
        )
        db_session.add(backup)
        db_session.commit()
        db_session.refresh(backup)
        
        # Mock S3 download (external service)
        mock_s3 = MagicMock()
        mock_s3_client.return_value = mock_s3
        mock_s3.generate_presigned_url.return_value = "https://s3.amazonaws.com/test-bucket/backup.sql?signed"
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/{backup.id}/download",
            "test_download_backup"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "download_url" in data
        assert "expires_at" in data
        assert data["download_url"].startswith("https://s3.amazonaws.com")

    def test_download_backup_not_completed(self, integration_client, db_session, timed_request):
        """Test downloading non-completed backup returns error."""
        # Create running backup in database
        backup = BackupRecord(
            id=uuid4(),
            backup_type="full",
            status="running",
            created_at=datetime.now(timezone.utc)
        )
        db_session.add(backup)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/{backup.id}/download",
            "test_download_not_completed"
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "not completed" in response.json()["detail"]

    @patch('amprenta_rag.api.routers.backup.export_project_data')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_create_project_export_success(self, mock_backup_engine, mock_export, 
                                         integration_client, db_session, test_program, timed_request):
        """Test creating project export with real database."""
        # Mock external services
        mock_task = MagicMock()
        mock_task.id = "export_task_123"
        mock_export.delay.return_value = mock_task
        
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/export/project",
            "test_create_project_export",
            json={
                "program_id": str(test_program.id),
                "entity_types": ["compounds", "datasets"],
                "format": "json"
            }
        )
        
        assert response.status_code == 202
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "export_id" in data
        assert data["message"] == "Export initiated"
        assert data["task_id"] == "export_task_123"
        
        # Verify export record was created in database
        export_id = data["export_id"]
        export_record = db_session.query(ProjectExport).filter_by(id=export_id).first()
        assert export_record is not None
        assert export_record.program_id == test_program.id
        assert export_record.status == "pending"

    def test_create_project_export_no_entities(self, integration_client, timed_request):
        """Test creating export with no entity types returns error."""
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/export/project",
            "test_export_no_entities",
            json={
                "program_id": str(uuid4()),
                "entity_types": [],
                "format": "json"
            }
        )
        
        assert response.status_code == 400
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "entity_types cannot be empty" in response.json()["detail"]

    def test_download_project_export_not_implemented(self, integration_client, timed_request):
        """Test project export download returns appropriate response."""
        export_id = uuid4()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{export_id}/download",
            "test_export_download_not_implemented"
        )
        
        # This endpoint may return 501 Not Implemented or 200 with implementation
        assert response.status_code in [200, 501]
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    @patch('amprenta_rag.api.routers.backup.export_project_data')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_create_returns_id(self, mock_backup_engine, mock_export_project,
                                    integration_client, db_session, test_program, timed_request):
        """Test export creation returns proper ID and database record."""
        # Mock external services
        mock_task = MagicMock()
        mock_task.id = "export_task_456"
        mock_export_project.delay.return_value = mock_task
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/export/project",
            "test_export_create_id",
            json={
                "program_id": str(test_program.id),
                "entity_types": ["compounds"],
                "format": "csv"
            }
        )
        
        assert response.status_code == 202
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        export_id = data["export_id"]
        
        # Verify database record exists
        export_record = db_session.query(ProjectExport).filter_by(id=export_id).first()
        assert export_record is not None
        assert export_record.format == "csv"
        assert len(export_record.entity_types) == 1
        assert "compounds" in export_record.entity_types

    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_download_success(self, mock_backup_engine, integration_client, 
                                   db_session, timed_request):
        """Test downloading completed export with mocked file system."""
        # Create completed export in database
        export = ProjectExport(
            id=uuid4(),
            program_id=uuid4(),
            entity_types=["compounds", "datasets"],
            format="json",
            status="completed",
            created_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            file_path="exports/test_export.zip",
            file_size=5000000
        )
        db_session.add(export)
        db_session.commit()
        db_session.refresh(export)
        
        # Mock backup engine (external file operations)
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.get_export_download_url.return_value = "https://download.example.com/export.zip"
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{export.id}/download",
            "test_export_download_success"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert "download_url" in data
        assert data["download_url"] == "https://download.example.com/export.zip"

    def test_export_download_expired_returns_410(self, integration_client, db_session, timed_request):
        """Test downloading expired export returns 410."""
        # Create expired export in database
        export = ProjectExport(
            id=uuid4(),
            program_id=uuid4(),
            entity_types=["compounds"],
            format="json",
            status="completed",
            created_at=datetime.now(timezone.utc) - timedelta(days=8),  # Expired
            completed_at=datetime.now(timezone.utc) - timedelta(days=8),
            file_path="exports/expired_export.zip",
            expires_at=datetime.now(timezone.utc) - timedelta(days=1)
        )
        db_session.add(export)
        db_session.commit()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{export.id}/download",
            "test_export_download_expired"
        )
        
        assert response.status_code == 410
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "expired" in response.json()["detail"].lower()

    def test_export_download_not_found_returns_404(self, integration_client, timed_request):
        """Test downloading non-existent export returns 404."""
        fake_export_id = uuid4()
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{fake_export_id}/download",
            "test_export_download_not_found"
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_cleanup_removes_expired(self, mock_backup_engine, integration_client, 
                                          db_session, timed_request):
        """Test cleanup endpoint removes expired exports from database."""
        # Create expired export in database
        expired_export = ProjectExport(
            id=uuid4(),
            program_id=uuid4(),
            entity_types=["compounds"],
            format="json",
            status="completed",
            created_at=datetime.now(timezone.utc) - timedelta(days=8),
            expires_at=datetime.now(timezone.utc) - timedelta(days=1)
        )
        
        # Create valid export in database
        valid_export = ProjectExport(
            id=uuid4(),
            program_id=uuid4(),
            entity_types=["datasets"],
            format="csv",
            status="completed",
            created_at=datetime.now(timezone.utc),
            expires_at=datetime.now(timezone.utc) + timedelta(days=6)
        )
        
        db_session.add_all([expired_export, valid_export])
        db_session.commit()
        
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.cleanup_expired_exports.return_value = 1
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/backup/export/cleanup",
            "test_export_cleanup"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        data = response.json()
        assert data["removed_count"] >= 1
        
        # Verify expired export was removed from database
        remaining_expired = db_session.query(ProjectExport).filter_by(id=expired_export.id).first()
        assert remaining_expired is None
        
        # Verify valid export remains
        remaining_valid = db_session.query(ProjectExport).filter_by(id=valid_export.id).first()
        assert remaining_valid is not None

    @patch('amprenta_rag.api.routers.backup.export_project_data')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_one_time_download(self, mock_backup_engine, mock_export_project,
                                    integration_client, db_session, test_program, timed_request):
        """Test one-time download functionality with database integration."""
        # Create export with one-time download flag
        export = ProjectExport(
            id=uuid4(),
            program_id=test_program.id,
            entity_types=["compounds"],
            format="json",
            status="completed",
            created_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            file_path="exports/one_time_export.zip",
            one_time_download=True
        )
        db_session.add(export)
        db_session.commit()
        db_session.refresh(export)
        
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.get_export_download_url.return_value = "https://download.example.com/one_time.zip"
        
        # First download should succeed
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{export.id}/download",
            "test_one_time_download_first"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        # Verify export was marked as downloaded in database
        db_session.refresh(export)
        assert export.downloaded_at is not None
        
        # Second download should fail
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/backup/export/{export.id}/download",
            "test_one_time_download_second"
        )
        
        assert response.status_code == 410
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        assert "already downloaded" in response.json()["detail"]
