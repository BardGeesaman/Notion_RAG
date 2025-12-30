"""Tests for sync Celery tasks."""

import os
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.sync import run_sync_job


class TestSyncTasks:
    """Test sync Celery tasks."""
    
    def setup_method(self):
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_success(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test successful sync task execution."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        
        # Mock job result
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.records_synced = 100
        mock_job.records_new = 50
        mock_job.records_updated = 30
        mock_job.conflicts_detected = 5
        mock_manager.run_sync.return_value = mock_job
        
        # Execute task
        result = run_sync_job(str(job_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["job_id"] == str(job_id)
        assert result["records_synced"] == 100
        assert result["records_new"] == 50
        assert result["records_updated"] == 30
        assert result["conflicts_detected"] == 5
        
        # Verify manager was created and adapters registered
        mock_manager_class.assert_called_once_with(mock_db_session)
        mock_manager.register_adapter.assert_any_call(mock_chembl)
        mock_manager.register_adapter.assert_any_call(mock_pubchem)
        mock_manager.run_sync.assert_called_once_with(job_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_job_not_found(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task when SyncJob doesn't exist."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager to raise exception for missing job
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        mock_manager.run_sync.side_effect = ValueError("Sync job not found")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            run_sync_job(str(job_id))
        
        # Verify correct error
        assert "Sync job not found" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_manager_error(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task when SyncManager raises exception."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager to raise exception
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        mock_manager.run_sync.side_effect = RuntimeError("Sync failed")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            run_sync_job(str(job_id))
        
        # Verify correct error
        assert "Sync failed" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_adapter_registration(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task verifies adapters are registered."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        
        # Mock job result
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.records_synced = 0
        mock_job.records_new = 0
        mock_job.records_updated = 0
        mock_job.conflicts_detected = 0
        mock_manager.run_sync.return_value = mock_job
        
        # Execute task
        result = run_sync_job(str(job_id))
        
        # Verify adapters were registered
        assert mock_manager.register_adapter.call_count == 2
        mock_manager.register_adapter.assert_any_call(mock_chembl)
        mock_manager.register_adapter.assert_any_call(mock_pubchem)
        
        # Verify adapter constructors were called correctly
        mock_chembl_adapter.assert_called_once()
        mock_pubchem_adapter.assert_called_once_with(mock_db_session)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_retry_on_failure(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task handles failure and re-raises for retry."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager to raise exception
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        mock_manager.run_sync.side_effect = Exception("Sync error")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            run_sync_job(str(job_id))
        
        # Verify correct error
        assert "Sync error" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_max_retries_exceeded(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task returns failure dict when max retries exceeded."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager to raise exception
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        mock_manager.run_sync.side_effect = Exception("Persistent sync failure")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            run_sync_job(str(job_id))
        
        # Verify correct error
        assert "Persistent sync failure" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_conflict_detection(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task returns conflicts_detected count."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        
        # Mock job result with conflicts
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.records_synced = 50
        mock_job.records_new = 20
        mock_job.records_updated = 25
        mock_job.conflicts_detected = 3
        mock_manager.run_sync.return_value = mock_job
        
        # Execute task
        result = run_sync_job(str(job_id))
        
        # Verify conflicts are reported
        assert result["status"] == "completed"
        assert result["conflicts_detected"] == 3
        assert result["records_synced"] == 50
        assert result["records_new"] == 20
        assert result["records_updated"] == 25
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.sync.manager.SyncManager')
    @patch('amprenta_rag.sync.adapters.chembl.ChEMBLAdapter')
    @patch('amprenta_rag.sync.adapters.pubchem.PubChemAdapter')
    def test_sync_db_update_on_failure(self, mock_pubchem_adapter, mock_chembl_adapter, mock_manager_class, mock_db_session):
        """Test sync task updates database status on failure."""
        job_id = uuid4()
        
        # Mock adapters
        mock_chembl = MagicMock()
        mock_pubchem = MagicMock()
        mock_chembl_adapter.return_value = mock_chembl
        mock_pubchem_adapter.return_value = mock_pubchem
        
        # Mock manager to raise exception
        mock_manager = MagicMock()
        mock_manager_class.return_value = mock_manager
        mock_manager.run_sync.side_effect = Exception("Sync processing error")
        
        # Mock job for database update
        mock_job = MagicMock()
        mock_job.status = "pending"
        mock_job.error_log = None
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = mock_job
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            run_sync_job(str(job_id))
        
        # Verify correct error
        assert "Sync processing error" in str(exc_info.value)
        
        # Verify job status was updated to failed
        assert mock_job.status == "failed"
        assert mock_job.completed_at is not None
        assert mock_job.error_log == "Sync processing error"
        mock_db.add.assert_called_once_with(mock_job)
        mock_db.commit.assert_called_once()
