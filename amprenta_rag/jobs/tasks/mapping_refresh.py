"""Celery tasks for ID mapping refresh and maintenance."""

import asyncio
from typing import Dict

from amprenta_rag.jobs.celery_app import celery_app
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@celery_app.task(bind=True, max_retries=3, default_retry_delay=300, queue='scheduled')
def refresh_uniprot_mappings_task(self) -> Dict[str, any]:
    """Weekly UniProt ID mapping refresh (bulk download).
    
    Downloads latest UniProt ID mappings and updates the database.
    Scheduled for Sunday 2am via Celery Beat.
    
    Returns:
        Dict with status and count of mappings updated
    """
    try:
        logger.info("Starting UniProt mappings refresh task")
        
        # Import here to avoid circular imports
        from amprenta_rag.services.id_mapping_service import refresh_uniprot_mappings
        
        # Run the async function
        count = asyncio.run(refresh_uniprot_mappings())
        
        logger.info(f"UniProt refresh completed: {count} mappings updated")
        return {"status": "success", "count": count}
        
    except Exception as exc:
        logger.exception("UniProt refresh failed")
        
        # Retry on failure with exponential backoff
        if self.request.retries < self.max_retries:
            countdown = 300 * (2 ** self.request.retries)  # 5, 10, 20 minutes
            logger.info(f"Retrying UniProt refresh (attempt {self.request.retries + 1}/{self.max_retries}) in {countdown}s")
            raise self.retry(exc=exc, countdown=countdown)
        
        # All retries exhausted
        logger.error(f"UniProt refresh failed after {self.max_retries} retries")
        return {"status": "failed", "error": str(exc), "retries": self.request.retries}


@celery_app.task(bind=True, queue='scheduled')
def cleanup_expired_mappings_task(self) -> Dict[str, any]:
    """Daily cleanup of expired KEGG cache entries.
    
    Removes mappings past their TTL. Scheduled for daily 3am.
    
    Returns:
        Dict with status and count of mappings removed
    """
    try:
        logger.info("Starting expired mappings cleanup task")
        
        # Import here to avoid circular imports
        from amprenta_rag.services.id_mapping_service import cleanup_expired_mappings
        
        count = cleanup_expired_mappings()
        
        logger.info(f"Expired mappings cleanup: {count} removed")
        return {"status": "success", "count": count}
        
    except Exception as exc:
        logger.exception("Cleanup failed")
        return {"status": "failed", "error": str(exc)}
