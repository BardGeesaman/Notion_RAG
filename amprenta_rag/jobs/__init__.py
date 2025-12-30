"""Job queue package for Celery-based task management."""

from .celery_app import celery_app

__all__ = ["celery_app"]
