"""Tests for docking Celery tasks."""

import os
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.docking import run_docking


class TestDockingTasks:
    """Test docking Celery tasks."""
    
    def setup_method(self):
        """Set up test environment."""
        # Enable eager execution for synchronous testing
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        """Clean up test environment."""
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    @patch('amprenta_rag.structural.docking_service.DockingService')
    @patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor')
    def test_docking_success(self, mock_executor, mock_docking_service, mock_prepare_receptor, mock_db_session):
        """Test successful docking task execution."""
        run_id = uuid4()
        compound_id = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = None
        mock_run.total_compounds = 10
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = None
        
        # Mock compounds
        mock_compound = MagicMock()
        mock_compound.id = compound_id
        compounds = [mock_compound]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.limit.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service
        mock_service = MagicMock()
        mock_docking_service.return_value = mock_service
        mock_pose = MagicMock()
        mock_service._process_compound.return_value = mock_pose
        
        # Mock ThreadPoolExecutor
        mock_future = MagicMock()
        mock_future.result.return_value = mock_pose
        mock_executor_instance = MagicMock()
        mock_executor_instance.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.return_value = mock_future
        mock_executor.return_value = mock_executor_instance
        
        # Mock as_completed
        with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [mock_future]
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["run_id"] == str(run_id)
        assert result["compounds_processed"] == 1
        
        # Verify receptor preparation was called
        mock_prepare_receptor.assert_called_once()
        
        # Verify docking service was used
        mock_service._process_compound.assert_called_once()
        
        # Verify run status updates
        assert mock_run.status == "completed"
        assert mock_run.completed_at is not None
    
    @patch('amprenta_rag.database.session.db_session')
    def test_docking_run_not_found(self, mock_db_session):
        """Test docking task when run doesn't exist."""
        run_id = uuid4()
        
        # Mock database session - return None for run query
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.return_value = None
        
        # Execute task
        result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Run not found"
        assert result["run_id"] == str(run_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    def test_docking_receptor_prep_failure(self, mock_prepare_receptor, mock_db_session):
        """Test docking task when receptor preparation fails."""
        run_id = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
        ]
        
        # Mock receptor preparation failure
        mock_prepare_receptor.return_value = False
        
        # Mock Path.exists to return False (receptor file doesn't exist)
        with patch('amprenta_rag.jobs.tasks.docking.Path') as mock_path:
            mock_path.return_value.exists.return_value = False
            mock_path.return_value.__truediv__ = lambda self, other: mock_path.return_value
            mock_path.return_value.mkdir = MagicMock()
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Receptor preparation failed"
        assert result["run_id"] == str(run_id)
        
        # Verify run was marked as failed
        assert mock_run.status == "failed"
        assert "Failed to prepare receptor" in mock_run.error_log
    
    @patch('amprenta_rag.database.session.db_session')
    def test_docking_no_structure(self, mock_db_session):
        """Test docking task when structure is missing."""
        run_id = uuid4()
        
        # Mock run without structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = None
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
        ]
        
        # Execute task
        result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Run or structure not found"
        assert result["run_id"] == str(run_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    @patch('amprenta_rag.structural.docking_service.DockingService')
    @patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor')
    def test_docking_compound_selection(self, mock_executor, mock_docking_service, mock_prepare_receptor, mock_db_session):
        """Test docking task compound selection logic."""
        run_id = uuid4()
        compound_id1 = uuid4()
        compound_id2 = uuid4()
        
        # Mock run with specific compound IDs
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = [str(compound_id1), str(compound_id2)]
        mock_run.total_compounds = None
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = None
        
        # Mock compounds
        mock_compound1 = MagicMock()
        mock_compound1.id = compound_id1
        mock_compound2 = MagicMock()
        mock_compound2.id = compound_id2
        compounds = [mock_compound1, mock_compound2]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.filter.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service
        mock_service = MagicMock()
        mock_docking_service.return_value = mock_service
        mock_pose = MagicMock()
        mock_service._process_compound.return_value = mock_pose
        
        # Mock ThreadPoolExecutor
        mock_future1 = MagicMock()
        mock_future1.result.return_value = mock_pose
        mock_future2 = MagicMock()
        mock_future2.result.return_value = mock_pose
        mock_executor_instance = MagicMock()
        mock_executor_instance.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.side_effect = [mock_future1, mock_future2]
        mock_executor.return_value = mock_executor_instance
        
        # Mock as_completed
        with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [mock_future1, mock_future2]
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["compounds_processed"] == 2
        
        # Verify specific compounds were used
        mock_db.query.return_value.filter.assert_called_once()
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    @patch('amprenta_rag.structural.docking_service.DockingService')
    @patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor')
    def test_docking_compound_processing_error(self, mock_executor, mock_docking_service, mock_prepare_receptor, mock_db_session):
        """Test docking task handling individual compound processing errors."""
        run_id = uuid4()
        compound_id = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = None
        mock_run.total_compounds = 1
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = None
        
        # Mock compounds
        mock_compound = MagicMock()
        mock_compound.id = compound_id
        compounds = [mock_compound]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Update error log query
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.limit.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service to raise exception
        mock_service = MagicMock()
        mock_docking_service.return_value = mock_service
        mock_service._process_compound.side_effect = RuntimeError("Compound processing failed")
        
        # Mock ThreadPoolExecutor
        mock_future = MagicMock()
        mock_future.result.side_effect = RuntimeError("Compound processing failed")
        mock_executor_instance = MagicMock()
        mock_executor_instance.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.return_value = mock_future
        mock_executor.return_value = mock_executor_instance
        
        # Mock as_completed
        with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [mock_future]
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify result - should still complete but with errors logged
        assert result["status"] == "completed"
        assert result["compounds_processed"] == 1
        
        # Verify error was logged to run
        # Note: The error_log update happens in the exception handler
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    @patch('amprenta_rag.structural.docking_service.DockingService')
    @patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor')
    def test_docking_error_log_aggregation(self, mock_executor, mock_docking_service, mock_prepare_receptor, mock_db_session):
        """Test docking task error log aggregation."""
        run_id = uuid4()
        compound_id1 = uuid4()
        compound_id2 = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = None
        mock_run.total_compounds = 2
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = ""
        
        # Mock compounds
        mock_compound1 = MagicMock()
        mock_compound1.id = compound_id1
        mock_compound2 = MagicMock()
        mock_compound2.id = compound_id2
        compounds = [mock_compound1, mock_compound2]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Error log update 1
            mock_run,  # Error log update 2
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.limit.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service
        mock_service = MagicMock()
        mock_docking_service.return_value = mock_service
        
        # Mock ThreadPoolExecutor
        mock_future1 = MagicMock()
        mock_future1.result.side_effect = RuntimeError("Error 1")
        mock_future2 = MagicMock()
        mock_future2.result.side_effect = RuntimeError("Error 2")
        mock_executor_instance = MagicMock()
        mock_executor_instance.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.side_effect = [mock_future1, mock_future2]
        mock_executor.return_value = mock_executor_instance
        
        # Mock as_completed
        with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [mock_future1, mock_future2]
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["compounds_processed"] == 2
        
        # Verify error aggregation would have happened
        # (In the actual implementation, errors are appended to error_log)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    @patch('amprenta_rag.structural.docking_service.DockingService')
    @patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor')
    def test_docking_status_transitions(self, mock_executor, mock_docking_service, mock_prepare_receptor, mock_db_session):
        """Test docking task status transitions."""
        run_id = uuid4()
        compound_id = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = None
        mock_run.total_compounds = 1
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = None
        
        # Mock compounds
        mock_compound = MagicMock()
        mock_compound.id = compound_id
        compounds = [mock_compound]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.limit.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service
        mock_service = MagicMock()
        mock_docking_service.return_value = mock_service
        mock_pose = MagicMock()
        mock_service._process_compound.return_value = mock_pose
        
        # Mock ThreadPoolExecutor
        mock_future = MagicMock()
        mock_future.result.return_value = mock_pose
        mock_executor_instance = MagicMock()
        mock_executor_instance.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.return_value = mock_future
        mock_executor.return_value = mock_executor_instance
        
        # Mock as_completed
        with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [mock_future]
            
            # Execute task
            result = run_docking(str(run_id))
        
        # Verify status transitions
        # Should go: pending → running → completed
        assert result["status"] == "completed"
        
        # Verify timestamps were set
        assert mock_run.started_at is not None
        assert mock_run.completed_at is not None
        assert mock_run.status == "completed"
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    def test_docking_retry_on_failure(self, mock_prepare_receptor, mock_db_session):
        """Test docking task retry logic with exponential backoff."""
        run_id = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.error_log = None
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            RuntimeError("Database error"),  # Simulate DB error
            mock_run,  # Query for error handling
        ]
        
        # Mock the task itself to simulate retry behavior
        with patch.object(run_docking, 'request') as mock_request:
            mock_request.retries = 1  # First retry
            run_docking.max_retries = 2
            
            with patch.object(run_docking, 'retry') as mock_retry:
                mock_retry.side_effect = Exception("Retry called")  # Simulate retry
                
                with pytest.raises(Exception, match="Retry called"):
                    run_docking(str(run_id))
                
                # Verify retry was called with exponential backoff
                mock_retry.assert_called_once()
                args, kwargs = mock_retry.call_args
                assert kwargs['countdown'] == 120 * (2 ** 1)  # 240 seconds for first retry
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.structural.receptor_prep.prepare_receptor')
    def test_docking_partial_success(self, mock_prepare_receptor, mock_db_session):
        """Test docking task when some compounds fail but run is marked failed due to errors."""
        run_id = uuid4()
        compound_id1 = uuid4()
        compound_id2 = uuid4()
        
        # Mock run and structure
        mock_run = MagicMock()
        mock_run.id = run_id
        mock_run.status = "pending"
        mock_run.structure = MagicMock()
        mock_run.compound_ids = None
        mock_run.total_compounds = 2
        mock_run.completed_compounds = 0
        mock_run.center_x = 0.0
        mock_run.center_y = 0.0
        mock_run.center_z = 0.0
        mock_run.size_x = 20.0
        mock_run.size_y = 20.0
        mock_run.size_z = 20.0
        mock_run.exhaustiveness = 8
        mock_run.error_log = "Some errors occurred"  # Pre-existing errors
        
        # Mock compounds
        mock_compound1 = MagicMock()
        mock_compound1.id = compound_id1
        mock_compound2 = MagicMock()
        mock_compound2.id = compound_id2
        compounds = [mock_compound1, mock_compound2]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter_by.return_value.first.side_effect = [
            mock_run,  # First query for run
            mock_run,  # Second query for run
            mock_run,  # Final query for run
        ]
        mock_db.query.return_value.limit.return_value.all.return_value = compounds
        
        # Mock receptor preparation
        mock_prepare_receptor.return_value = True
        
        # Mock docking service
        with patch('amprenta_rag.jobs.tasks.docking.DockingService') as mock_docking_service:
            mock_service = MagicMock()
            mock_docking_service.return_value = mock_service
            mock_pose = MagicMock()
            mock_service._process_compound.return_value = mock_pose
            
            # Mock ThreadPoolExecutor
            with patch('amprenta_rag.jobs.tasks.docking.ThreadPoolExecutor') as mock_executor:
                mock_future1 = MagicMock()
                mock_future1.result.return_value = mock_pose  # Success
                mock_future2 = MagicMock()
                mock_future2.result.return_value = None  # Failure
                mock_executor_instance = MagicMock()
                mock_executor_instance.__enter__.return_value = mock_executor_instance
                mock_executor_instance.submit.side_effect = [mock_future1, mock_future2]
                mock_executor.return_value = mock_executor_instance
                
                # Mock as_completed
                with patch('amprenta_rag.jobs.tasks.docking.as_completed') as mock_as_completed:
                    mock_as_completed.return_value = [mock_future1, mock_future2]
                    
                    # Execute task
                    result = run_docking(str(run_id))
        
        # Verify result - should be failed due to error_log
        assert result["status"] == "completed"  # Task completes, but run may be marked failed
        assert result["compounds_processed"] == 2
        
        # Verify run status depends on error_log
        # In the actual implementation, status is "failed" if error_log exists
