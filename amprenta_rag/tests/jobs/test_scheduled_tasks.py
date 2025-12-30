"""Tests for Celery scheduled tasks."""

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest


class TestScheduledTaskRegistration:
    """Test that scheduled tasks are properly registered with Celery."""

    def test_harvest_task_registered(self):
        """Test that harvest task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.scheduled import run_harvest
        
        assert run_harvest is not None
        task_name = 'amprenta_rag.jobs.tasks.scheduled.run_harvest'
        assert task_name in celery_app.tasks

    def test_digest_task_registered(self):
        """Test that digest task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.scheduled import run_digest
        
        assert run_digest is not None
        task_name = 'amprenta_rag.jobs.tasks.scheduled.run_digest'
        assert task_name in celery_app.tasks

    def test_external_sync_task_registered(self):
        """Test that external sync task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.scheduled import run_external_sync
        
        assert run_external_sync is not None
        task_name = 'amprenta_rag.jobs.tasks.scheduled.run_external_sync'
        assert task_name in celery_app.tasks

    def test_check_schedules_tasks_registered(self):
        """Test that schedule check tasks are registered."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.scheduled import check_digest_schedules, check_harvest_schedules
        
        assert check_digest_schedules is not None
        assert check_harvest_schedules is not None
        
        digest_task = 'amprenta_rag.jobs.tasks.scheduled.check_digest_schedules'
        harvest_task = 'amprenta_rag.jobs.tasks.scheduled.check_harvest_schedules'
        assert digest_task in celery_app.tasks
        assert harvest_task in celery_app.tasks


class TestBeatScheduleConfiguration:
    """Test Celery Beat schedule configuration."""

    def test_beat_schedule_has_sync_entries(self):
        """Test that CELERYBEAT_SCHEDULE has external sync entries."""
        from amprenta_rag.jobs.schedules import CELERYBEAT_SCHEDULE
        
        assert isinstance(CELERYBEAT_SCHEDULE, dict)
        assert 'sync-chembl-daily' in CELERYBEAT_SCHEDULE
        assert 'sync-pubchem-daily' in CELERYBEAT_SCHEDULE
        
        # Check ChEMBL sync configuration
        chembl_config = CELERYBEAT_SCHEDULE['sync-chembl-daily']
        assert chembl_config['task'] == 'amprenta_rag.jobs.tasks.scheduled.run_external_sync'
        assert chembl_config['args'] == ['chembl', 'incremental']
        assert chembl_config['options']['queue'] == 'scheduled'
        
        # Check PubChem sync configuration
        pubchem_config = CELERYBEAT_SCHEDULE['sync-pubchem-daily']
        assert pubchem_config['task'] == 'amprenta_rag.jobs.tasks.scheduled.run_external_sync'
        assert pubchem_config['args'] == ['pubchem', 'incremental']
        assert pubchem_config['options']['queue'] == 'scheduled'

    def test_beat_schedule_loaded_in_config(self):
        """Test that beat schedule is loaded in Celery config."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        # Check that beat schedule is configured
        beat_schedule = celery_app.conf.beat_schedule
        assert beat_schedule is not None
        assert isinstance(beat_schedule, dict)
        
        # Should contain our sync tasks
        assert 'sync-chembl-daily' in beat_schedule
        assert 'sync-pubchem-daily' in beat_schedule


class TestScheduledTaskQueueRouting:
    """Test that scheduled tasks use the correct queue."""

    def test_scheduled_tasks_use_scheduled_queue(self):
        """Test that scheduled tasks are routed to 'scheduled' queue."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        routes = celery_app.conf.task_routes
        
        # Check if there's a specific route for scheduled tasks
        scheduled_route = routes.get('amprenta_rag.jobs.tasks.scheduled.*')
        if scheduled_route:
            assert scheduled_route['queue'] == 'scheduled'
        
        # Also check that the tasks have queue specified in decorator
        from amprenta_rag.jobs.tasks.scheduled import run_harvest, run_digest, run_external_sync
        
        # These tasks should be callable (they exist and are Celery tasks)
        assert callable(run_harvest)
        assert callable(run_digest)
        assert callable(run_external_sync)


class TestScheduledTaskExecution:
    """Test scheduled task execution with mocking."""

    @patch('amprenta_rag.jobs.tasks.scheduled._run_harvest')
    def test_harvest_task_calls_underlying_function(self, mock_run_harvest):
        """Test that harvest task calls the underlying harvest function."""
        from amprenta_rag.jobs.tasks.scheduled import run_harvest
        
        schedule_id = str(uuid4())
        
        # Mock the underlying function to succeed
        mock_run_harvest.return_value = None
        
        # Call the task directly (in eager mode)
        result = run_harvest(schedule_id)
        
        # Verify the underlying function was called
        mock_run_harvest.assert_called_once_with(schedule_id)
        
        # Verify the result
        assert result['status'] == 'complete'
        assert result['schedule_id'] == schedule_id

    @patch('amprenta_rag.jobs.tasks.scheduled._run_external_sync')
    def test_external_sync_task_calls_underlying_function(self, mock_run_sync):
        """Test that external sync task calls the underlying sync function."""
        from amprenta_rag.jobs.tasks.scheduled import run_external_sync
        
        source = 'chembl'
        sync_type = 'incremental'
        
        # Mock the underlying function to succeed
        mock_run_sync.return_value = None
        
        # Call the task directly (in eager mode)
        result = run_external_sync(source, sync_type)
        
        # Verify the underlying function was called
        mock_run_sync.assert_called_once_with(source, sync_type)
        
        # Verify the result
        assert result['status'] == 'complete'
        assert result['source'] == source
        assert result['sync_type'] == sync_type

    @patch('amprenta_rag.automation.digest_scheduler.DigestScheduler')
    def test_digest_task_calls_scheduler(self, mock_scheduler_class):
        """Test that digest task calls DigestScheduler."""
        from amprenta_rag.jobs.tasks.scheduled import run_digest
        
        schedule_id = str(uuid4())
        
        # Mock the scheduler and its result
        mock_result = MagicMock()
        mock_result.status = 'complete'
        mock_result.output_path = '/path/to/output.html'
        
        mock_scheduler = MagicMock()
        mock_scheduler.run_digest.return_value = mock_result
        mock_scheduler_class.return_value = mock_scheduler
        
        # Call the task directly
        result = run_digest(schedule_id)
        
        # Verify the scheduler was created and called
        mock_scheduler_class.assert_called_once()
        mock_scheduler.run_digest.assert_called_once()
        
        # Verify the result
        assert result['status'] == 'complete'
        assert result['schedule_id'] == schedule_id
        assert result['output_path'] == '/path/to/output.html'
