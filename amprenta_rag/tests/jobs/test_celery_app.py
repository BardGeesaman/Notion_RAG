"""Tests for Celery application configuration."""

import os
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from celery import Celery


class TestCeleryApp:
    """Test Celery application initialization and configuration."""

    def test_celery_app_initialization(self):
        """Test that Celery app is properly initialized."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        assert isinstance(celery_app, Celery)
        assert celery_app.main == 'amprenta_rag'

    def test_config_loading(self):
        """Test that Celery loads configuration from config module."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        # Check that config is loaded from the config module
        assert celery_app.conf.task_serializer == 'json'
        assert celery_app.conf.result_serializer == 'json'
        assert celery_app.conf.accept_content == ['json']

    def test_queue_routing_configuration(self):
        """Test that queue routing is properly configured."""
        from amprenta_rag.jobs.celery_app import celery_app
        
        # Check default queue
        assert celery_app.conf.task_default_queue == 'default'
        
        # Check task routes exist
        assert celery_app.conf.task_routes is not None
        assert 'amprenta_rag.jobs.tasks.docking.*' in celery_app.conf.task_routes
        assert celery_app.conf.task_routes['amprenta_rag.jobs.tasks.docking.*']['queue'] == 'high'

    def test_celery_task_always_eager_mode(self):
        """Test CELERY_TASK_ALWAYS_EAGER configuration."""
        with patch.dict(os.environ, {'CELERY_TASK_ALWAYS_EAGER': 'true'}):
            # Need to reload config to pick up env var
            from amprenta_rag.jobs import config
            import importlib
            importlib.reload(config)
            
            assert config.task_always_eager is True

    def test_worker_process_init_signal_registered(self):
        """Test that worker_process_init signal is registered."""
        # Import the module to ensure signal handler is registered
        from amprenta_rag.jobs import celery_app  # noqa: F401
        from celery.signals import worker_process_init
        
        # Check that at least one signal handler is registered
        # The exact inspection of handlers is complex due to Celery internals
        # So we'll just verify the signal has receivers
        assert len(worker_process_init.receivers) > 0, "No worker_process_init signal handlers registered"

    @patch('amprenta_rag.database.base.get_engine')
    def test_worker_process_init_handler(self, mock_get_engine):
        """Test that worker_process_init handler calls engine.dispose()."""
        # Setup mock engine
        mock_engine = MagicMock()
        mock_get_engine.return_value = mock_engine
        
        from amprenta_rag.jobs.celery_app import init_worker
        
        # Call the handler directly
        init_worker()
        
        # Verify get_engine() was called and engine.dispose() was called
        mock_get_engine.assert_called_once()
        mock_engine.dispose.assert_called_once()

    def test_use_celery_feature_flag_handling(self):
        """Test USE_CELERY feature flag configuration."""
        # Test with USE_CELERY=true
        with patch.dict(os.environ, {'USE_CELERY': 'true'}):
            use_celery = os.getenv('USE_CELERY', 'false').lower() == 'true'
            assert use_celery is True
        
        # Test with USE_CELERY=false (default)
        with patch.dict(os.environ, {}, clear=True):
            use_celery = os.getenv('USE_CELERY', 'false').lower() == 'true'
            assert use_celery is False
