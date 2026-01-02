"""Celery tasks for system maintenance and cleanup."""

from celery import shared_task
from amprenta_rag.services.lifecycle import cleanup_orphans
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@shared_task(name="cleanup_orphans_task")
def cleanup_orphans_task(dry_run: bool = False) -> dict:
    """
    Scheduled task to clean up orphaned entities.
    
    Run weekly: celery -A amprenta_rag.jobs beat
    """
    logger.info(f"Starting orphan cleanup (dry_run={dry_run})")
    results = cleanup_orphans(dry_run=dry_run)
    logger.info(f"Orphan cleanup complete: {results}")
    return results
