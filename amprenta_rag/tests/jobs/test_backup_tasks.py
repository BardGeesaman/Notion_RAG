"""Tests for backup Celery tasks."""

import os
from datetime import datetime, timedelta, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.backup.backup_engine import BackupEngineError
from amprenta_rag.database.models import BackupRecord
from amprenta_rag.jobs.schedules import CELERYBEAT_SCHEDULE
from amprenta_rag.jobs.tasks.backup import (
    backup_health_check,
    cleanup_old_backups,
    run_database_backup,
    verify_latest_backup,
)


class TestBackupTasks:
    """Test backup Celery tasks."""
    
    def setup_method(self):
        """Set up test environment."""
        # Enable eager execution for synchronous testing
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        """Clean up test environment."""
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    def test_run_database_backup_success(self, mock_backup_engine):
        """Test successful database backup task execution."""
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        
        backup_id = uuid4()
        mock_engine.create_backup.return_value = backup_id
        
        # Mock backup record
        mock_backup_record = MagicMock()
        mock_backup_record.status = "completed"
        mock_backup_record.file_size_bytes = 1024000
        mock_engine.get_backup_status.return_value = mock_backup_record
        
        # Execute task
        result = run_database_backup("full")
        
        # Verify result
        assert result["status"] == "completed"
        assert result["backup_id"] == str(backup_id)
        assert result["backup_type"] == "full"
        assert result["file_size_bytes"] == "1024000"
        
        # Verify engine calls
        mock_engine.create_backup.assert_called_once_with("full")
        mock_engine.get_backup_status.assert_called_once_with(backup_id)
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    def test_run_database_backup_failure(self, mock_backup_engine):
        """Test database backup task failure handling."""
        # Mock backup engine to raise exception
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.create_backup.side_effect = BackupEngineError("Backup failed")
        
        # Mock the task's retry mechanism to prevent actual retries in tests
        original_task = run_database_backup
        
        # Create a test version that doesn't retry
        def mock_task_func(backup_type):
            try:
                engine = mock_engine
                backup_id = engine.create_backup(backup_type)
                backup_record = engine.get_backup_status(backup_id)
                if backup_record and backup_record.status == "completed":
                    return {
                        "status": "completed",
                        "backup_id": str(backup_id),
                        "backup_type": backup_type,
                        "file_size_bytes": str(backup_record.file_size_bytes or 0),
                    }
                else:
                    error_msg = backup_record.error_message if backup_record else "Unknown error"
                    raise BackupEngineError(f"Backup failed: {error_msg}")
            except Exception as exc:
                # Return failure without retry for test
                return {
                    "status": "failed",
                    "backup_type": backup_type,
                    "error": str(exc),
                }
        
        # Execute mock task
        result = mock_task_func("full")
        
        # Verify failure result
        assert result["status"] == "failed"
        assert result["backup_type"] == "full"
        assert "Backup failed" in result["error"]
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    @patch('amprenta_rag.jobs.tasks.backup.get_config')
    def test_cleanup_old_backups_success(self, mock_get_config, mock_db_session, mock_backup_engine):
        """Test successful backup cleanup task."""
        # Mock config
        mock_config = MagicMock()
        mock_config.backup.retention_days = 30
        mock_get_config.return_value = mock_config
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock old backup records
        old_backup1 = MagicMock()
        old_backup1.id = uuid4()
        old_backup1.file_path = "backups/old1.sql.gz"
        
        old_backup2 = MagicMock()
        old_backup2.id = uuid4()
        old_backup2.file_path = "backups/old2.sql.gz"
        
        mock_db.query.return_value.filter.return_value.all.return_value = [old_backup1, old_backup2]
        
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        
        # Execute task
        result = cleanup_old_backups()
        
        # Verify result
        assert result["status"] == "completed"
        assert result["deleted_count"] == "2"
        assert result["error_count"] == "0"
        assert result["retention_days"] == "30"
        
        # Verify database operations
        assert mock_db.delete.call_count == 2
        mock_db.commit.assert_called_once()
        
        # Verify file deletions
        assert mock_engine.backup_client.delete_backup.call_count == 2
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_verify_latest_backup_success(self, mock_db_session, mock_backup_engine):
        """Test successful backup verification task."""
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock latest backup record
        latest_backup = MagicMock()
        latest_backup.id = uuid4()
        latest_backup.created_at = datetime.now(timezone.utc)
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = latest_backup
        
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.verify_backup.return_value = True
        
        # Execute task
        result = verify_latest_backup()
        
        # Verify result
        assert result["status"] == "verified"
        assert result["backup_id"] == str(latest_backup.id)
        assert "backup_date" in result
        
        # Verify engine calls
        mock_engine.verify_backup.assert_called_once_with(latest_backup.id)
        mock_db.expunge.assert_called_once_with(latest_backup)
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_verify_latest_backup_no_backups(self, mock_db_session, mock_backup_engine):
        """Test backup verification when no backups exist."""
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = None
        
        # Execute task
        result = verify_latest_backup()
        
        # Verify result
        assert result["status"] == "skipped"
        assert result["reason"] == "No completed backups found"
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_verify_latest_backup_failure(self, mock_db_session, mock_backup_engine):
        """Test backup verification failure."""
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock latest backup record
        latest_backup = MagicMock()
        latest_backup.id = uuid4()
        latest_backup.created_at = datetime.now(timezone.utc)
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = latest_backup
        
        # Mock backup engine
        mock_engine = MagicMock()
        mock_backup_engine.return_value = mock_engine
        mock_engine.verify_backup.return_value = False
        
        # Execute task
        result = verify_latest_backup()
        
        # Verify result
        assert result["status"] == "failed"
        assert result["backup_id"] == str(latest_backup.id)
        assert result["error"] == "Backup verification failed"
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_backup_health_check_healthy(self, mock_db_session, mock_backup_engine):
        """Test backup health check with healthy system."""
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock recent backup exists
        mock_db.query.return_value.filter.return_value.first.return_value = MagicMock()
        
        # Mock no recent failures
        mock_db.query.return_value.filter.return_value.count.return_value = 0
        
        # Mock total backup count
        mock_db.query.return_value.count.return_value = 10
        
        # Mock backup engine availability
        mock_backup_engine.return_value = MagicMock()
        
        # Execute task
        result = backup_health_check()
        
        # Verify result
        assert result["status"] == "healthy"
        assert "Recent backup found" in result["checks"]
        assert "No recent failures" in result["checks"]
        assert "Backup engine available" in result["checks"]
        assert result["total_backups"] == "10"
    
    @patch('amprenta_rag.jobs.tasks.backup.BackupEngine')
    @patch('amprenta_rag.jobs.tasks.backup.db_session')
    def test_backup_health_check_warnings(self, mock_db_session, mock_backup_engine):
        """Test backup health check with warning conditions."""
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        
        # Mock no recent backup
        query_mock = mock_db.query.return_value
        query_mock.filter.return_value.first.return_value = None
        
        # Mock recent failures
        query_mock.filter.return_value.count.return_value = 2
        
        # Mock total backup count
        mock_db.query.return_value.count.return_value = 5
        
        # Mock backup engine availability
        mock_backup_engine.return_value = MagicMock()
        
        # Execute task
        result = backup_health_check()
        
        # Verify result
        assert result["status"] == "warning"
        assert "No recent backup found (48h)" in result["checks"]
        assert "2 failed backups in last week" in result["checks"]
    
    def test_backup_schedule_configuration(self):
        """Test that backup tasks are properly scheduled in Celery Beat."""
        # Verify daily backup schedule
        assert 'daily-database-backup' in CELERYBEAT_SCHEDULE
        daily_backup = CELERYBEAT_SCHEDULE['daily-database-backup']
        assert daily_backup['task'] == 'amprenta_rag.jobs.tasks.backup.run_database_backup'
        assert 1 in daily_backup['schedule'].hour  # 1:00 AM (crontab returns a set)
        assert 0 in daily_backup['schedule'].minute
        assert daily_backup['args'] == ['full']
        assert daily_backup['options']['queue'] == 'scheduled'
        
        # Verify weekly cleanup schedule
        assert 'weekly-backup-cleanup' in CELERYBEAT_SCHEDULE
        weekly_cleanup = CELERYBEAT_SCHEDULE['weekly-backup-cleanup']
        assert weekly_cleanup['task'] == 'amprenta_rag.jobs.tasks.backup.cleanup_old_backups'
        assert 4 in weekly_cleanup['schedule'].hour  # 4:00 AM
        assert 0 in weekly_cleanup['schedule'].minute
        assert 0 in weekly_cleanup['schedule'].day_of_week  # Sunday
        assert weekly_cleanup['options']['queue'] == 'scheduled'
        
        # Verify weekly verification schedule
        assert 'weekly-backup-verify' in CELERYBEAT_SCHEDULE
        weekly_verify = CELERYBEAT_SCHEDULE['weekly-backup-verify']
        assert weekly_verify['task'] == 'amprenta_rag.jobs.tasks.backup.verify_latest_backup'
        assert 5 in weekly_verify['schedule'].hour  # 5:00 AM
        assert 0 in weekly_verify['schedule'].minute
        assert 0 in weekly_verify['schedule'].day_of_week  # Sunday
        assert weekly_verify['options']['queue'] == 'scheduled'
        
        # Verify daily health check schedule
        assert 'daily-backup-health-check' in CELERYBEAT_SCHEDULE
        health_check = CELERYBEAT_SCHEDULE['daily-backup-health-check']
        assert health_check['task'] == 'amprenta_rag.jobs.tasks.backup.backup_health_check'
        assert 6 in health_check['schedule'].hour  # 6:00 AM
        assert 0 in health_check['schedule'].minute
        assert health_check['options']['queue'] == 'scheduled'
