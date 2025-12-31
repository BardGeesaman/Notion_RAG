"""Tests for backup API endpoints."""

import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.api.main import app
from amprenta_rag.database.models import BackupRecord, ProjectExport
from amprenta_rag.models.auth import User


class TestBackupAPI:
    """Test backup API endpoints."""
    
    def setup_method(self):
        """Set up test environment."""
        # Enable eager execution for synchronous testing
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
        
        # Create test client
        self.client = TestClient(app)
        
        # Mock user for authentication
        self.test_user_id = uuid4()
        self.test_user = User(
            id=self.test_user_id,
            username="testuser",
            email="test@example.com",
            password_hash="fake_hash"
        )
    
    def teardown_method(self):
        """Clean up test environment."""
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
        # Clear any remaining dependency overrides
        app.dependency_overrides.clear()
    
    def _get_auth_headers(self) -> dict:
        """Get authentication headers for test user."""
        return {"X-User-Id": str(self.test_user_id)}
    
    @patch('amprenta_rag.api.routers.backup.run_database_backup')
    @patch('amprenta_rag.api.dependencies.get_database_session')
    def test_trigger_manual_backup_success(self, mock_get_db, mock_backup_task):
        """Test successful manual backup trigger."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        mock_get_db.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = self.test_user
        
        # Mock Celery task
        mock_task = MagicMock()
        mock_task.id = "task_12345"
        mock_backup_task.delay.return_value = mock_task
        
        # Override dependencies
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            response = self.client.post(
                "/api/v1/backup/database",
                json={"backup_type": "full"},
                headers=self._get_auth_headers(),
            )
            
            # Verify response
            assert response.status_code == 202
            data = response.json()
            assert data["message"] == "Full backup initiated"
            assert data["task_id"] == "task_12345"
            assert data["backup_type"] == "full"
            
            # Verify task was called
            mock_backup_task.delay.assert_called_once_with("full")
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_trigger_manual_backup_invalid_type(self):
        """Test manual backup with invalid backup type."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        
        try:
            # Make request with invalid type
            response = self.client.post(
                "/api/v1/backup/database",
                json={"backup_type": "invalid"},
                headers=self._get_auth_headers(),
            )
            
            # Verify error response
            assert response.status_code == 400
            assert "backup_type must be" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_trigger_manual_backup_unauthorized(self):
        """Test manual backup without authentication."""
        response = self.client.post(
            "/api/v1/backup/database",
            json={"backup_type": "full"},
        )
        
        # Verify unauthorized response
        assert response.status_code == 401
    
    def test_get_backup_history_success(self):
        """Test successful backup history retrieval."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        
        # Mock backup records
        backup1 = BackupRecord(
            id=uuid4(),
            backup_type="full",
            status="completed",
            file_path="backups/backup1.sql.gz",
            file_size_bytes=1024000,
            checksum_sha256="abc123",
            started_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            created_at=datetime.now(timezone.utc),
        )
        
        backup2 = BackupRecord(
            id=uuid4(),
            backup_type="incremental",
            status="running",
            started_at=datetime.now(timezone.utc),
            created_at=datetime.now(timezone.utc),
        )
        
        mock_backups = [backup1, backup2]
        
        # Mock query chain
        mock_query = mock_db.query.return_value
        mock_query.count.return_value = 2
        mock_query.order_by.return_value.offset.return_value.limit.return_value.all.return_value = mock_backups
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            response = self.client.get(
                "/api/v1/backup/history?page=1&per_page=10",
                headers=self._get_auth_headers(),
            )
            
            # Verify response
            assert response.status_code == 200
            data = response.json()
            assert data["total"] == 2
            assert data["page"] == 1
            assert data["per_page"] == 10
            assert len(data["items"]) == 2
            
            # Verify first backup
            item1 = data["items"][0]
            assert item1["backup_type"] == "full"
            assert item1["status"] == "completed"
            assert item1["file_size_bytes"] == 1024000
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_get_backup_details_success(self):
        """Test successful backup details retrieval."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        
        # Mock backup record
        backup_id = uuid4()
        backup = BackupRecord(
            id=backup_id,
            backup_type="full",
            status="completed",
            file_path="backups/backup.sql.gz",
            file_size_bytes=2048000,
            checksum_sha256="def456",
            started_at=datetime.now(timezone.utc),
            completed_at=datetime.now(timezone.utc),
            created_at=datetime.now(timezone.utc),
        )
        
        mock_db.query.return_value.filter.return_value.first.return_value = backup
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            response = self.client.get(
                f"/api/v1/backup/{backup_id}",
                headers=self._get_auth_headers(),
            )
            
            # Verify response
            assert response.status_code == 200
            data = response.json()
            assert data["id"] == str(backup_id)
            assert data["backup_type"] == "full"
            assert data["status"] == "completed"
            assert data["file_size_bytes"] == 2048000
            assert data["checksum_sha256"] == "def456"
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_get_backup_details_not_found(self):
        """Test backup details for non-existent backup."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            backup_id = uuid4()
            response = self.client.get(
                f"/api/v1/backup/{backup_id}",
                headers=self._get_auth_headers(),
            )
            
            # Verify error response
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_download_backup_success(self):
        """Test backup download endpoint validation."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        
        # Mock backup record (not found to test validation)
        backup_id = uuid4()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request - should get 404 for not found
            response = self.client.get(
                f"/api/v1/backup/{backup_id}/download",
                headers=self._get_auth_headers(),
            )
            
            # Should get a 404 error for backup not found
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_download_backup_not_completed(self):
        """Test download backup that is not completed."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        
        # Mock backup record (not completed)
        backup_id = uuid4()
        backup = BackupRecord(
            id=backup_id,
            backup_type="full",
            status="running",
            file_path="backups/backup.sql.gz",
            created_at=datetime.now(timezone.utc),
        )
        
        mock_db.query.return_value.filter.return_value.first.return_value = backup
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            response = self.client.get(
                f"/api/v1/backup/{backup_id}/download",
                headers=self._get_auth_headers(),
            )
            
            # Verify error response
            assert response.status_code == 400
            assert "not completed" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    @patch('amprenta_rag.api.routers.backup.export_project')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_create_project_export_success(self, mock_backup_engine, mock_export):
        """Test successful project export creation."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        # Mock database session
        mock_db = MagicMock()
        
        # Mock export function
        mock_export_data = b"test export data"
        mock_export.return_value = mock_export_data
        
        # Mock backup engine
        mock_client = MagicMock()
        mock_backup_engine.return_value.backup_client = mock_client
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        app.dependency_overrides[get_database_session] = lambda: mock_db
        
        try:
            # Make request
            program_ids = [uuid4(), uuid4()]
            experiment_ids = [uuid4()]
            
            response = self.client.post(
                "/api/v1/backup/export",
                json={
                    "program_ids": [str(pid) for pid in program_ids],
                    "experiment_ids": [str(eid) for eid in experiment_ids],
                    "include_related": True,
                },
                headers=self._get_auth_headers(),
            )
            
            # Verify response (now 201 instead of 202)
            assert response.status_code == 201
            data = response.json()
            assert data["message"] == "Project export created successfully"
            assert data["export_size_bytes"] == len(mock_export_data)
            assert "programs: 2" in data["entities_summary"]
            assert "experiments: 1" in data["entities_summary"]
            assert "compounds: 0" in data["entities_summary"]
            assert "export_id" in data
            assert "expires_at" in data
            
            # Verify export function was called
            mock_export.assert_called_once_with(
                program_ids=program_ids,
                experiment_ids=experiment_ids,
                compound_ids=None,
                db=mock_db,
                include_related=True,
            )
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_create_project_export_no_entities(self):
        """Test project export with no entities specified."""
        # Override authentication dependency
        def mock_get_current_user():
            return self.test_user
        
        app.dependency_overrides[get_current_user] = mock_get_current_user
        
        try:
            # Make request with no entities
            response = self.client.post(
                "/api/v1/backup/export",
                json={"include_related": True},
                headers=self._get_auth_headers(),
            )
            
            # Verify error response
            assert response.status_code == 400
            assert "At least one of" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    def test_download_project_export_now_implemented(self):
        """Test project export download is now implemented (returns 404 for non-existent)."""
        try:
            # Mock database session with no export found
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            export_id = uuid4()
            
            # Make request
            response = self.client.get(
                f"/api/v1/backup/export/{export_id}",
                headers=self._get_auth_headers(),
            )
            
            # Verify 404 response (no longer 501)
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]
            
        finally:
            # Clear dependency overrides
            app.dependency_overrides.clear()
    
    @patch('amprenta_rag.api.routers.backup.export_project')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_create_returns_id(self, mock_backup_engine, mock_export_project):
        """Test that POST /export creates export and returns ID."""
        try:
            # Mock database session
            mock_db = MagicMock()
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            # Mock export_project to return ZIP data
            mock_export_data = b"fake zip data"
            mock_export_project.return_value = mock_export_data
            
            # Mock backup engine and client
            mock_client = MagicMock()
            mock_backup_engine.return_value.backup_client = mock_client
            
            # Make request
            response = self.client.post(
                "/api/v1/backup/export",
                json={
                    "program_ids": [str(uuid4())],
                    "include_related": True
                },
                headers=self._get_auth_headers(),
            )
            
            # Verify response
            assert response.status_code == 201
            data = response.json()
            assert "export_id" in data
            assert data["message"] == "Project export created successfully"
            assert data["export_size_bytes"] == len(mock_export_data)
            assert "expires_at" in data
            
            # Verify database operations
            mock_db.add.assert_called_once()
            mock_db.commit.assert_called_once()
            
        finally:
            app.dependency_overrides.clear()
    
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_download_success(self, mock_backup_engine):
        """Test successful export download."""
        try:
            # Mock database session and ProjectExport
            mock_db = MagicMock()
            export_id = uuid4()
            
            from datetime import datetime, timedelta
            from amprenta_rag.database.models import ProjectExport
            
            mock_export = ProjectExport(
                id=export_id,
                file_path="exports/test.zip",
                file_size_bytes=1024,
                entities_summary="test",
                created_by=self.test_user_id,
                expires_at=datetime.utcnow() + timedelta(hours=1),  # Not expired
                created_at=datetime.utcnow()
            )
            
            mock_db.query.return_value.filter.return_value.first.return_value = mock_export
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            # Mock backup engine and file operations
            mock_client = MagicMock()
            mock_backup_engine.return_value.backup_client = mock_client
            
            # Mock temporary file creation and file operations
            with patch('amprenta_rag.api.routers.backup.tempfile.NamedTemporaryFile') as mock_temp:
                with patch('amprenta_rag.api.routers.backup.Path') as mock_path:
                    # Create a real temporary file for testing
                    with tempfile.NamedTemporaryFile(delete=False) as real_temp:
                        real_temp.write(b"test zip content")
                        real_temp.flush()
                        temp_path = real_temp.name
                    
                    # Mock the temporary file creation
                    mock_temp.return_value.__enter__.return_value.name = temp_path
                    mock_path.return_value.exists.return_value = True
                    
                    try:
                        # Make request
                        response = self.client.get(
                            f"/api/v1/backup/export/{export_id}",
                            headers=self._get_auth_headers(),
                        )
                    finally:
                        # Clean up the real temp file
                        import os
                        if os.path.exists(temp_path):
                            os.unlink(temp_path)
            
            # Verify response
            assert response.status_code == 200
            assert response.headers["content-type"] == "application/zip"
            
            # Verify download and cleanup
            mock_client.download_file.assert_called_once()
            mock_db.delete.assert_called_once_with(mock_export)
            mock_db.commit.assert_called()
            
        finally:
            app.dependency_overrides.clear()
    
    def test_export_download_expired_returns_410(self):
        """Test that downloading expired export returns 410 GONE."""
        try:
            # Mock database session with expired export
            mock_db = MagicMock()
            export_id = uuid4()
            
            from datetime import datetime, timedelta
            from amprenta_rag.database.models import ProjectExport
            
            mock_export = ProjectExport(
                id=export_id,
                file_path="exports/test.zip",
                file_size_bytes=1024,
                entities_summary="test",
                created_by=self.test_user_id,
                expires_at=datetime.utcnow() - timedelta(hours=1),  # Expired
                created_at=datetime.utcnow()
            )
            
            mock_db.query.return_value.filter.return_value.first.return_value = mock_export
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            # Make request
            response = self.client.get(
                f"/api/v1/backup/export/{export_id}",
                headers=self._get_auth_headers(),
            )
            
            # Verify 410 GONE response
            assert response.status_code == 410
            assert "expired" in response.json()["detail"]
            
            # Verify expired export was deleted
            mock_db.delete.assert_called_once_with(mock_export)
            mock_db.commit.assert_called_once()
            
        finally:
            app.dependency_overrides.clear()
    
    def test_export_download_not_found_returns_404(self):
        """Test that downloading non-existent export returns 404."""
        try:
            # Mock database session with no export found
            mock_db = MagicMock()
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            export_id = uuid4()
            
            # Make request
            response = self.client.get(
                f"/api/v1/backup/export/{export_id}",
                headers=self._get_auth_headers(),
            )
            
            # Verify 404 response
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]
            
        finally:
            app.dependency_overrides.clear()
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_export_cleanup_removes_expired(self, mock_db_session, mock_backup_engine):
        """Test that cleanup task removes expired exports."""
        from datetime import datetime, timedelta
        from amprenta_rag.database.models import ProjectExport
        from amprenta_rag.jobs.tasks.backup import cleanup_expired_exports
        
        # Mock database session context manager
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Create mock expired exports
        expired_export1 = MagicMock(spec=ProjectExport)
        expired_export1.id = uuid4()
        expired_export1.file_path = "exports/expired1.zip"
        
        expired_export2 = MagicMock(spec=ProjectExport)
        expired_export2.id = uuid4()
        expired_export2.file_path = "exports/expired2.zip"
        
        mock_db.query.return_value.filter.return_value.all.return_value = [
            expired_export1, expired_export2
        ]
        
        # Mock backup engine
        mock_client = MagicMock()
        mock_backup_engine.return_value.backup_client = mock_client
        
        # Run cleanup task
        result = cleanup_expired_exports()
        
        # Verify results
        assert result["status"] == "completed"
        assert result["deleted_count"] == "2"
        assert result["error_count"] == "0"
        
        # Verify file deletions
        assert mock_client.delete_backup.call_count == 2
        mock_client.delete_backup.assert_any_call("exports/expired1.zip")
        mock_client.delete_backup.assert_any_call("exports/expired2.zip")
        
        # Verify database deletions
        assert mock_db.delete.call_count == 2
        mock_db.commit.assert_called_once()
    
    @patch('amprenta_rag.api.routers.backup.export_project')
    @patch('amprenta_rag.api.routers.backup.BackupEngine')
    def test_export_one_time_download(self, mock_backup_engine, mock_export_project):
        """Test that export can only be downloaded once."""
        try:
            # Mock database session
            mock_db = MagicMock()
            app.dependency_overrides[get_database_session] = lambda: mock_db
            app.dependency_overrides[get_current_user] = lambda: self.test_user
            
            # First, create an export
            mock_export_data = b"fake zip data"
            mock_export_project.return_value = mock_export_data
            
            mock_client = MagicMock()
            mock_backup_engine.return_value.backup_client = mock_client
            
            # Create export
            response = self.client.post(
                "/api/v1/backup/export",
                json={
                    "program_ids": [str(uuid4())],
                    "include_related": True
                },
                headers=self._get_auth_headers(),
            )
            
            assert response.status_code == 201
            export_id = response.json()["export_id"]
            
            # Reset mock to simulate first download
            mock_db.reset_mock()
            
            from datetime import datetime, timedelta
            from amprenta_rag.database.models import ProjectExport
            
            mock_export = ProjectExport(
                id=export_id,
                file_path="exports/test.zip",
                file_size_bytes=1024,
                entities_summary="test",
                created_by=self.test_user_id,
                expires_at=datetime.utcnow() + timedelta(hours=1),
                created_at=datetime.utcnow()
            )
            
            # Mock first download - export exists
            mock_db.query.return_value.filter.return_value.first.return_value = mock_export
            
            # Mock temporary file creation and file operations
            with patch('amprenta_rag.api.routers.backup.tempfile.NamedTemporaryFile') as mock_temp:
                with patch('amprenta_rag.api.routers.backup.Path') as mock_path:
                    # Create a real temporary file for testing
                    with tempfile.NamedTemporaryFile(delete=False) as real_temp:
                        real_temp.write(b"test zip content")
                        real_temp.flush()
                        temp_path = real_temp.name
                    
                    # Mock the temporary file creation
                    mock_temp.return_value.__enter__.return_value.name = temp_path
                    mock_path.return_value.exists.return_value = True
                    
                    try:
                        # First download should succeed
                        response = self.client.get(
                            f"/api/v1/backup/export/{export_id}",
                            headers=self._get_auth_headers(),
                        )
                    finally:
                        # Clean up the real temp file
                        import os
                        if os.path.exists(temp_path):
                            os.unlink(temp_path)
            
            assert response.status_code == 200
            
            # Verify export was deleted after download
            mock_db.delete.assert_called_once_with(mock_export)
            
            # Reset mock for second attempt
            mock_db.reset_mock()
            
            # Mock second download - export no longer exists
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            # Second download should return 404
            response = self.client.get(
                f"/api/v1/backup/export/{export_id}",
                headers=self._get_auth_headers(),
            )
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"]
            
        finally:
            app.dependency_overrides.clear()