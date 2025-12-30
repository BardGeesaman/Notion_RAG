"""Tests for Celery configuration validation."""

import importlib
import os
from unittest.mock import patch

import pytest


class TestBrokerConfiguration:
    """Test broker and backend configuration."""

    def test_broker_url_default(self):
        """Test default broker URL."""
        with patch.dict(os.environ, {}, clear=True):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.broker_url == 'redis://localhost:6379/0'

    def test_broker_url_from_env(self):
        """Test broker URL from environment variable."""
        with patch.dict(os.environ, {'CELERY_BROKER_URL': 'redis://custom:6380/1'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.broker_url == 'redis://custom:6380/1'

    def test_result_backend_default(self):
        """Test default result backend."""
        with patch.dict(os.environ, {}, clear=True):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.result_backend == 'redis://localhost:6379/0'

    def test_result_backend_from_env(self):
        """Test result backend from environment variable."""
        with patch.dict(os.environ, {'CELERY_RESULT_BACKEND': 'redis://results:6381/2'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.result_backend == 'redis://results:6381/2'


class TestSerializationSettings:
    """Test serialization configuration."""

    def test_task_serializer_is_json(self):
        """Test that task serializer is JSON."""
        from amprenta_rag.jobs import config
        
        assert config.task_serializer == 'json'

    def test_result_serializer_is_json(self):
        """Test that result serializer is JSON."""
        from amprenta_rag.jobs import config
        
        assert config.result_serializer == 'json'

    def test_accept_content_json_only(self):
        """Test that only JSON content is accepted."""
        from amprenta_rag.jobs import config
        
        assert config.accept_content == ['json']
        assert 'pickle' not in config.accept_content  # Security: no pickle


class TestResultSettings:
    """Test result expiration and storage settings."""

    def test_result_expires_default(self):
        """Test default result expiration (24 hours)."""
        with patch.dict(os.environ, {}, clear=True):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.result_expires == 86400  # 24 hours in seconds

    def test_result_expires_from_env(self):
        """Test result expiration from environment variable."""
        with patch.dict(os.environ, {'CELERY_RESULT_EXPIRES': '3600'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.result_expires == 3600


class TestWorkerSettings:
    """Test worker configuration."""

    def test_worker_prefetch_multiplier(self):
        """Test worker prefetch multiplier is 1 for long-running tasks."""
        from amprenta_rag.jobs import config
        
        assert config.worker_prefetch_multiplier == 1

    def test_worker_concurrency_default(self):
        """Test default worker concurrency."""
        with patch.dict(os.environ, {}, clear=True):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.worker_concurrency == 4

    def test_worker_concurrency_from_env(self):
        """Test worker concurrency from environment variable."""
        with patch.dict(os.environ, {'CELERY_WORKER_CONCURRENCY': '8'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.worker_concurrency == 8


class TestQueueDefinitions:
    """Test queue configuration."""

    def test_default_queue_exists(self):
        """Test that default queue is defined."""
        from amprenta_rag.jobs import config
        
        queue_names = [q.name for q in config.task_queues]
        assert 'default' in queue_names

    def test_high_priority_queue_exists(self):
        """Test that high priority queue is defined."""
        from amprenta_rag.jobs import config
        
        queue_names = [q.name for q in config.task_queues]
        assert 'high' in queue_names

    def test_low_priority_queue_exists(self):
        """Test that low priority queue is defined."""
        from amprenta_rag.jobs import config
        
        queue_names = [q.name for q in config.task_queues]
        assert 'low' in queue_names

    def test_scheduled_queue_exists(self):
        """Test that scheduled queue is defined."""
        from amprenta_rag.jobs import config
        
        queue_names = [q.name for q in config.task_queues]
        assert 'scheduled' in queue_names

    def test_queue_count(self):
        """Test that all 4 queues are defined."""
        from amprenta_rag.jobs import config
        
        assert len(config.task_queues) == 4


class TestTaskRouting:
    """Test task routing configuration."""

    def test_genomics_route_to_default(self):
        """Test genomics tasks route to default queue."""
        from amprenta_rag.jobs import config
        
        route = config.task_routes.get('amprenta_rag.jobs.tasks.genomics.*')
        assert route is not None
        assert route['queue'] == 'default'

    def test_docking_route_to_high(self):
        """Test docking tasks route to high priority queue."""
        from amprenta_rag.jobs import config
        
        route = config.task_routes.get('amprenta_rag.jobs.tasks.docking.*')
        assert route is not None
        assert route['queue'] == 'high'

    def test_extraction_route_to_high(self):
        """Test extraction tasks route to high priority queue."""
        from amprenta_rag.jobs import config
        
        route = config.task_routes.get('amprenta_rag.jobs.tasks.extraction.*')
        assert route is not None
        assert route['queue'] == 'high'

    def test_sync_route_to_low(self):
        """Test sync tasks route to low priority queue."""
        from amprenta_rag.jobs import config
        
        route = config.task_routes.get('amprenta_rag.jobs.tasks.sync.*')
        assert route is not None
        assert route['queue'] == 'low'

    def test_single_cell_route_to_default(self):
        """Test single cell tasks route to default queue."""
        from amprenta_rag.jobs import config
        
        route = config.task_routes.get('amprenta_rag.jobs.tasks.single_cell.*')
        assert route is not None
        assert route['queue'] == 'default'

    def test_default_queue_setting(self):
        """Test default queue is 'default'."""
        from amprenta_rag.jobs import config
        
        assert config.task_default_queue == 'default'


class TestBeatScheduleConfiguration:
    """Test Celery Beat schedule configuration."""

    def test_beat_schedule_imported(self):
        """Test that beat schedule is imported from schedules module."""
        from amprenta_rag.jobs import config
        
        assert config.beat_schedule is not None
        assert isinstance(config.beat_schedule, dict)

    def test_beat_scheduler_type(self):
        """Test that persistent scheduler is used."""
        from amprenta_rag.jobs import config
        
        assert config.beat_scheduler == 'celery.beat:PersistentScheduler'


class TestEagerModeConfiguration:
    """Test task eager mode configuration for testing."""

    def test_eager_mode_default_false(self):
        """Test that eager mode is disabled by default."""
        with patch.dict(os.environ, {}, clear=True):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.task_always_eager is False

    def test_eager_mode_enabled_from_env(self):
        """Test that eager mode can be enabled via environment."""
        with patch.dict(os.environ, {'CELERY_TASK_ALWAYS_EAGER': 'true'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.task_always_eager is True

    def test_eager_mode_case_insensitive(self):
        """Test that eager mode env var is case insensitive."""
        with patch.dict(os.environ, {'CELERY_TASK_ALWAYS_EAGER': 'TRUE'}):
            from amprenta_rag.jobs import config
            importlib.reload(config)
            
            assert config.task_always_eager is True

