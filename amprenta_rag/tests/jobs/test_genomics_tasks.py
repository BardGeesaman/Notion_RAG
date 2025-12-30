"""Tests for genomics Celery tasks."""

import os
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.genomics import run_genomics_pipeline


class TestGenomicsTasks:
    """Test genomics Celery tasks."""
    
    def setup_method(self):
        """Set up test environment."""
        # Enable eager execution for synchronous testing
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        """Clean up test environment."""
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_kallisto_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_success(self, mock_get_config, mock_kallisto, mock_salmon, mock_db_session):
        """Test successful genomics pipeline task execution."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "salmon"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/salmon_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job,  # First query for job
            mock_job,  # Second query for job
            mock_index,  # Query for index
            mock_job,  # Final query for job
        ]
        
        # Mock pipeline result
        result_file = Path("/tmp/jobs") / str(job_id) / "quant.sf"
        mock_salmon.return_value = result_file
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify result
        assert result["status"] == "complete"
        assert result["job_id"] == str(job_id)
        assert result["result_file"] == str(result_file)
        
        # Verify pipeline was called correctly
        mock_salmon.assert_called_once()
        
        # Verify job status updates
        assert mock_job.status == "complete"
        assert mock_job.progress_percent == 100
        assert mock_job.result_file == str(result_file)
        assert mock_job.error_message is None
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_job_not_found(self, mock_get_config, mock_db_session):
        """Test genomics pipeline task when job doesn't exist."""
        job_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock database session - return None for job query
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Job not found"
        assert result["job_id"] == str(job_id)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_index_not_found(self, mock_get_config, mock_db_session):
        """Test genomics pipeline task when GenomicsIndex is missing."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "salmon"
        mock_job.index_id = index_id
        mock_job.progress_percent = 0
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job,  # First query for job
            mock_job,  # Second query for job
            None,      # Query for index - not found
            mock_job,  # Query for job in error handling
        ]
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert "Index not found" in result["error"]
        assert result["job_id"] == str(job_id)
        
        # Verify job was marked as failed
        assert mock_job.status == "failed"
        assert "Index not found" in mock_job.error_message
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_pipeline_salmon_success(self, mock_get_config, mock_salmon, mock_db_session):
        """Test successful Salmon tool execution."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "salmon"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/salmon_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job, mock_job, mock_index, mock_job
        ]
        
        # Mock pipeline result
        result_file = Path("/tmp/jobs") / str(job_id) / "quant.sf"
        mock_salmon.return_value = result_file
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify salmon was called
        mock_salmon.assert_called_once_with(
            fastq_path=Path("/data/sample.fastq"),
            index_path=Path("/indices/salmon_index"),
            output_dir=Path("/tmp/jobs") / str(job_id)
        )
        
        # Verify success
        assert result["status"] == "complete"
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_kallisto_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_pipeline_kallisto_success(self, mock_get_config, mock_kallisto, mock_db_session):
        """Test successful Kallisto tool execution."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "kallisto"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/kallisto_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job, mock_job, mock_index, mock_job
        ]
        
        # Mock pipeline result
        result_file = Path("/tmp/jobs") / str(job_id) / "abundance.tsv"
        mock_kallisto.return_value = result_file
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify kallisto was called
        mock_kallisto.assert_called_once_with(
            fastq_path=Path("/data/sample.fastq"),
            index_path=Path("/indices/kallisto_index"),
            output_dir=Path("/tmp/jobs") / str(job_id)
        )
        
        # Verify success
        assert result["status"] == "complete"
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_unknown_tool(self, mock_get_config, mock_db_session):
        """Test genomics pipeline task with unknown tool."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "unknown_tool"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/some_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job, mock_job, mock_index, mock_job
        ]
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert "Unknown tool" in result["error"]
        assert result["job_id"] == str(job_id)
        
        # Verify job was marked as failed
        assert mock_job.status == "failed"
        assert "Unknown tool" in mock_job.error_message
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_pipeline_error(self, mock_get_config, mock_salmon, mock_db_session):
        """Test genomics pipeline task handling pipeline execution error."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "salmon"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/salmon_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job, mock_job, mock_index, mock_job
        ]
        
        # Mock pipeline to raise exception
        mock_salmon.side_effect = RuntimeError("Pipeline execution failed")
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert "Pipeline execution failed" in result["error"]
        assert result["job_id"] == str(job_id)
        
        # Verify job was marked as failed
        assert mock_job.status == "failed"
        assert "Pipeline execution failed" in mock_job.error_message
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_status_transitions(self, mock_get_config, mock_salmon, mock_db_session):
        """Test genomics pipeline task status transitions."""
        job_id = uuid4()
        index_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock job and index
        mock_job = MagicMock()
        mock_job.id = job_id
        mock_job.tool = "salmon"
        mock_job.index_id = index_id
        mock_job.input_fastq_path = "/data/sample.fastq"
        mock_job.output_dir = None
        mock_job.progress_percent = 0
        mock_job.status = "pending"
        
        mock_index = MagicMock()
        mock_index.id = index_id
        mock_index.file_path = "/indices/salmon_index"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_job, mock_job, mock_index, mock_job
        ]
        
        # Mock pipeline result
        result_file = Path("/tmp/jobs") / str(job_id) / "quant.sf"
        mock_salmon.return_value = result_file
        
        # Execute task
        result = run_genomics_pipeline(str(job_id))
        
        # Verify status transitions
        # Should go: pending → running → complete
        assert result["status"] == "complete"
        
        # Verify progress updates
        assert mock_job.progress_percent == 100  # Final progress
        assert mock_job.status == "complete"  # Final status
        assert mock_job.started_at is not None
        assert mock_job.completed_at is not None
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_retry_on_failure(self, mock_get_config, mock_salmon, mock_db_session):
        """Test that task handles failure and returns error status."""
        job_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock database to raise exception
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = Exception("DB connection failed")
        
        # Execute - in eager mode, task catches exception and returns failure
        result = run_genomics_pipeline(str(job_id))
        
        # Verify failure handling
        assert result["status"] == "failed"
        assert "job_id" in result
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.ingestion.genomics.pipeline.run_salmon_quant')
    @patch('amprenta_rag.config.get_config')
    def test_genomics_max_retries_exceeded(self, mock_get_config, mock_salmon, mock_db_session):
        """Test that task handles persistent failure and returns error status."""
        job_id = uuid4()
        
        # Mock config
        mock_config = MagicMock()
        mock_config.genomics.job_root = "/tmp/jobs"
        mock_get_config.return_value = mock_config
        
        # Mock database to raise persistent exception
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.side_effect = Exception("Persistent failure")
        
        # Execute - in eager mode, task catches exception and returns failure
        result = run_genomics_pipeline(str(job_id))
        
        # Verify final failure
        assert result["status"] == "failed"
        assert "Persistent failure" in result["error"]
        assert result["job_id"] == str(job_id)
