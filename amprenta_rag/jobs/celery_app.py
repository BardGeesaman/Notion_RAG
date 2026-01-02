"""Celery application configuration."""

from celery import Celery
from celery.signals import worker_process_init, worker_ready

celery_app = Celery('amprenta_rag')
celery_app.config_from_object('amprenta_rag.jobs.config')


@worker_process_init.connect
def init_worker(**kwargs):
    """Initialize worker process for SQLAlchemy fork safety."""
    from amprenta_rag.database.base import get_engine
    engine = get_engine()
    engine.dispose()  # Force new connections in worker


@worker_ready.connect
def verify_worker_secrets(sender, **kwargs):
    """Verify worker has required secrets on startup."""
    from amprenta_rag.utils.config_check import validate_required_secrets
    import logging
    
    logger = logging.getLogger(__name__)
    secrets_valid, missing = validate_required_secrets()
    if not secrets_valid:
        logger.critical(f"Celery worker missing secrets: {missing}")
        # Don't crash worker, but log critical error
        # In production, monitoring should alert on this
