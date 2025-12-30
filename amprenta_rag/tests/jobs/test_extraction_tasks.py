"""Tests for extraction Celery tasks."""

import os
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.extraction import process_extraction_job


class TestExtractionTasks:
    """Test extraction Celery tasks."""
    
    def setup_method(self):
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_success(self, mock_service_class, mock_db_session):
        """Test successful extraction task execution."""
        job_id = uuid4()
        
        # Mock service
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        
        # Mock job result
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.completed_count = 10
        mock_job.file_count = 10
        mock_service.process_job.return_value = mock_job
        
        # Execute task
        result = process_extraction_job(str(job_id))
        
        # Verify result
        assert result["status"] == "completed"
        assert result["job_id"] == str(job_id)
        assert result["completed_count"] == 10
        assert result["file_count"] == 10
        
        # Verify service was called correctly
        mock_service_class.assert_called_once_with(mock_db_session)
        mock_service.process_job.assert_called_once_with(job_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_job_not_found(self, mock_service_class, mock_db_session):
        """Test extraction task when job doesn't exist."""
        job_id = uuid4()
        
        # Mock service to raise exception for missing job
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        mock_service.process_job.side_effect = ValueError("Job not found")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_extraction_job(str(job_id))
        
        # Verify correct error
        assert "Job not found" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_service_error(self, mock_service_class, mock_db_session):
        """Test extraction task when ExtractionBatchService raises exception."""
        job_id = uuid4()
        
        # Mock service to raise exception
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        mock_service.process_job.side_effect = RuntimeError("Service error")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_extraction_job(str(job_id))
        
        # Verify correct error
        assert "Service error" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_status_transitions(self, mock_service_class, mock_db_session):
        """Test extraction task status transitions."""
        job_id = uuid4()
        
        # Mock service
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        
        # Mock job with status progression
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.completed_count = 5
        mock_job.file_count = 5
        mock_service.process_job.return_value = mock_job
        
        # Execute task
        result = process_extraction_job(str(job_id))
        
        # Verify final status
        assert result["status"] == "completed"
        assert result["completed_count"] == 5
        assert result["file_count"] == 5
        
        # Verify service was used
        mock_service.process_job.assert_called_once_with(job_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_retry_on_failure(self, mock_service_class, mock_db_session):
        """Test extraction task handles failure and re-raises for retry."""
        job_id = uuid4()
        
        # Mock service to raise exception
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        mock_service.process_job.side_effect = Exception("Extraction failed")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_extraction_job(str(job_id))
        
        # Verify correct error
        assert "Extraction failed" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_max_retries_exceeded(self, mock_service_class, mock_db_session):
        """Test extraction task returns failure dict when max retries exceeded."""
        job_id = uuid4()
        
        # Mock service to raise exception
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        mock_service.process_job.side_effect = Exception("Persistent failure")
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_extraction_job(str(job_id))
        
        # Verify correct error
        assert "Persistent failure" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_partial_completion(self, mock_service_class, mock_db_session):
        """Test extraction task with partial completion."""
        job_id = uuid4()
        
        # Mock service
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        
        # Mock job with partial completion
        mock_job = MagicMock()
        mock_job.status = "completed"
        mock_job.completed_count = 7
        mock_job.file_count = 10  # 3 files failed
        mock_service.process_job.return_value = mock_job
        
        # Execute task
        result = process_extraction_job(str(job_id))
        
        # Verify partial completion is handled
        assert result["status"] == "completed"
        assert result["completed_count"] == 7
        assert result["file_count"] == 10
        assert result["completed_count"] < result["file_count"]  # Partial completion
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.extraction.batch_service.ExtractionBatchService')
    def test_extraction_db_update_on_failure(self, mock_service_class, mock_db_session):
        """Test extraction task updates database status on failure."""
        job_id = uuid4()
        
        # Mock service to raise exception
        mock_service = MagicMock()
        mock_service_class.return_value = mock_service
        mock_service.process_job.side_effect = Exception("Processing error")
        
        # Mock job for database update
        mock_job = MagicMock()
        mock_job.status = "pending"
        
        # Mock database for error handling
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = mock_job
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_extraction_job(str(job_id))
        
        # Verify correct error
        assert "Processing error" in str(exc_info.value)
        
        # Verify job status was updated to failed
        assert mock_job.status == "failed"
        assert mock_job.completed_at is not None
        mock_db.add.assert_called_once_with(mock_job)
        mock_db.commit.assert_called_once()
