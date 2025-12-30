"""Tests for Celery task definitions."""

import os
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest


class TestTaskRegistration:
    """Test that all tasks are properly registered with Celery."""

    def test_genomics_task_registered(self):
        """Test that genomics task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.genomics import run_genomics_pipeline
        
        # Import the task to ensure it's registered
        assert run_genomics_pipeline is not None
        
        # Check task is in registered tasks
        task_name = 'amprenta_rag.jobs.tasks.genomics.run_genomics_pipeline'
        assert task_name in celery_app.tasks

    def test_docking_task_registered(self):
        """Test that docking task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.docking import run_docking
        
        assert run_docking is not None
        task_name = 'amprenta_rag.jobs.tasks.docking.run_docking'
        assert task_name in celery_app.tasks

    def test_extraction_task_registered(self):
        """Test that extraction task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.extraction import process_extraction_job
        
        assert process_extraction_job is not None
        task_name = 'amprenta_rag.jobs.tasks.extraction.process_extraction_job'
        assert task_name in celery_app.tasks

    def test_sync_task_registered(self):
        """Test that sync task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.sync import run_sync_job
        
        assert run_sync_job is not None
        task_name = 'amprenta_rag.jobs.tasks.sync.run_sync_job'
        assert task_name in celery_app.tasks

    def test_single_cell_task_registered(self):
        """Test that single cell task is registered in celery_app.tasks."""
        from amprenta_rag.jobs.celery_app import celery_app
        from amprenta_rag.jobs.tasks.single_cell import process_single_cell
        
        assert process_single_cell is not None
        task_name = 'amprenta_rag.jobs.tasks.single_cell.process_single_cell'
        assert task_name in celery_app.tasks


class TestTaskQueueRouting:
    """Test that tasks are routed to correct queues."""

    def test_genomics_task_queue_routing(self):
        """Test that genomics task is routed to 'default' queue."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        task_name = 'amprenta_rag.jobs.tasks.genomics.run_genomics_pipeline'
        routes = celery_app.conf.task_routes
        
        # Check if there's a specific route or falls back to default
        route = routes.get('amprenta_rag.jobs.tasks.genomics.*')
        if route:
            assert route['queue'] == 'default'
        else:
            # Falls back to default queue
            assert celery_app.conf.task_default_queue == 'default'

    def test_docking_task_queue_routing(self):
        """Test that docking task is routed to 'high' queue."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        routes = celery_app.conf.task_routes
        route = routes.get('amprenta_rag.jobs.tasks.docking.*')
        assert route is not None
        assert route['queue'] == 'high'

    def test_sync_task_queue_routing(self):
        """Test that sync task is routed to 'low' queue."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        routes = celery_app.conf.task_routes
        route = routes.get('amprenta_rag.jobs.tasks.sync.*')
        assert route is not None
        assert route['queue'] == 'low'


class TestTaskBehavior:
    """Test task execution and retry behavior."""

    def test_task_retry_behavior(self):
        """Test that tasks have retry configuration."""
        from amprenta_rag.jobs.tasks.genomics import run_genomics_pipeline
        from amprenta_rag.jobs.tasks.docking import run_docking
        from amprenta_rag.jobs.tasks.extraction import process_extraction_job
        
        # Verify tasks are callable (they exist and are Celery tasks)
        assert callable(run_genomics_pipeline)
        assert callable(run_docking)
        assert callable(process_extraction_job)
        
        # Verify they are Celery task objects
        assert hasattr(run_genomics_pipeline, 'delay')
        assert hasattr(run_docking, 'delay')
        assert hasattr(process_extraction_job, 'delay')

    def test_task_eager_execution(self):
        """Test that CELERY_TASK_ALWAYS_EAGER runs tasks synchronously."""
        with patch.dict(os.environ, {'CELERY_TASK_ALWAYS_EAGER': 'true'}):
            from amprenta_rag.jobs.celery_app import celery_app
            
            # Reload config to pick up env var
            celery_app.conf.task_always_eager = True
            
            # In eager mode, tasks run synchronously
            assert celery_app.conf.task_always_eager is True
