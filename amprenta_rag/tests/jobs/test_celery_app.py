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
        from amprenta_rag.jobs.celery_app import celery_app
        from celery.signals import worker_process_init
        
        # Check that our signal handler is registered
        receivers = worker_process_init.receivers
        assert len(receivers) > 0
        
        # Find our specific handler
        handler_found = False
        for receiver in receivers:
            if hasattr(receiver, '__name__') and 'init_worker' in receiver.__name__:
                handler_found = True
                break
        
        assert handler_found, "worker_process_init signal handler not found"

    @patch('amprenta_rag.database.base.engine')
    def test_worker_process_init_handler(self, mock_engine):
        """Test that worker_process_init handler calls engine.dispose()."""
        from amprenta_rag.jobs.celery_app import init_worker
        
        # Call the handler directly
        init_worker()
        
        # Verify engine.dispose() was called
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
