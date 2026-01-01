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


@celery_app.task(bind=True, name="cleanup_expired_mappings")
def cleanup_expired_mappings_task(self) -> dict:
    """Daily task to remove expired ID mappings."""
    try:
        logger.info("Starting expired mappings cleanup task")
        
        # Import here to avoid circular imports
        from amprenta_rag.services.id_mapping_service import cleanup_expired_mappings
        
        deleted_count = cleanup_expired_mappings()
        
        logger.info(f"Expired mappings cleanup: {deleted_count} removed")
        return {
            "status": "success",
            "deleted_count": deleted_count,
        }
        
    except Exception as exc:
        logger.exception("Cleanup failed")
        return {"status": "failed", "error": str(exc)}


@celery_app.task(bind=True, name="refresh_kegg_cache")
def refresh_kegg_cache_task(self) -> dict:
    """Daily task to refresh expiring KEGG mappings."""
    import asyncio
    
    try:
        logger.info("Starting KEGG cache refresh task")
        
        # Import here to avoid circular imports
        from amprenta_rag.sync.adapters.kegg_refresh import KEGGRefreshAdapter
        from amprenta_rag.services.id_mapping_service import (
            log_refresh_start,
            log_refresh_complete,
        )
        
        log_entry = log_refresh_start("kegg_refresh")
        
        try:
            adapter = KEGGRefreshAdapter()
            refreshed_count = 0
            
            async def _run():
                nonlocal refreshed_count
                async for record in adapter.fetch_records():
                    if record.get("refreshed"):
                        refreshed_count += 1
            
            asyncio.run(_run())
            
            log_refresh_complete(str(log_entry.id), refreshed_count)
            
            logger.info(f"KEGG cache refresh completed: {refreshed_count} mappings refreshed")
            return {
                "status": "success",
                "refreshed_count": refreshed_count,
                "log_id": str(log_entry.id),
            }
            
        except Exception as e:
            log_refresh_complete(str(log_entry.id), 0, error=str(e))
            raise
            
    except Exception as exc:
        logger.exception("KEGG cache refresh failed")
        return {"status": "failed", "error": str(exc)}


@celery_app.task(bind=True, name="auto_resolve_sync_conflicts")
def auto_resolve_sync_conflicts_task(self, source: str = None, limit: int = 100) -> dict:
    """Task to auto-resolve pending sync conflicts."""
    try:
        logger.info(f"Starting auto-resolve sync conflicts task (source={source}, limit={limit})")
        
        # Import here to avoid circular imports
        from amprenta_rag.sync.conflict_resolver import auto_resolve_pending_conflicts
        
        resolved_count = auto_resolve_pending_conflicts(source=source, limit=limit)
        
        logger.info(f"Auto-resolve conflicts completed: {resolved_count} conflicts resolved")
        return {
            "status": "success",
            "resolved_count": resolved_count,
            "source": source,
        }
        
    except Exception as exc:
        logger.exception("Auto-resolve conflicts failed")
        return {"status": "failed", "error": str(exc)}
