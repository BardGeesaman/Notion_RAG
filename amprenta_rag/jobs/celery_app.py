"""Celery application configuration."""

from celery import Celery
from celery.signals import worker_process_init

celery_app = Celery('amprenta_rag')
celery_app.config_from_object('amprenta_rag.jobs.config')


@worker_process_init.connect
def init_worker(**kwargs):
    """Initialize worker process for SQLAlchemy fork safety."""
    from amprenta_rag.database.base import get_engine
    engine = get_engine()
    engine.dispose()  # Force new connections in worker
