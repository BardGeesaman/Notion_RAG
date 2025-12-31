"""Admin API endpoints for system management."""

from fastapi import APIRouter, Depends, HTTPException
from typing import Dict, Any

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.admin.cache_manager import get_cache_stats, clear_cache, clear_all_caches, get_cache_summary
from amprenta_rag.utils.health import get_extended_system_info, get_celery_queue_stats, get_connection_status
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)
router = APIRouter(prefix="/admin", tags=["admin"])


def _require_admin(current_user) -> None:
    """Helper to ensure current user is admin."""
    if not current_user or current_user.role != "admin":
        raise HTTPException(status_code=403, detail="Admin access required")


@router.get("/caches")
def list_caches(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """List all caches with detailed statistics (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s requested cache stats", current_user.id)
    return get_cache_stats()


@router.get("/caches/summary")
def get_caches_summary(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Get high-level cache usage summary (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s requested cache summary", current_user.id)
    return get_cache_summary()


@router.post("/caches/{cache_name}/clear")
def clear_specific_cache(cache_name: str, current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Clear a specific cache (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s clearing cache: %s", current_user.id, cache_name)
    success = clear_cache(cache_name)
    
    if not success:
        raise HTTPException(
            status_code=404, 
            detail=f"Cache '{cache_name}' not found or failed to clear"
        )
    
    return {
        "status": "cleared",
        "cache": cache_name,
        "message": f"Successfully cleared {cache_name}"
    }


@router.post("/caches/clear-all")
def clear_all(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Clear all caches (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s clearing all caches", current_user.id)
    results = clear_all_caches()
    
    # Count successes and failures
    successful = sum(1 for success in results.values() if success)
    failed = len(results) - successful
    
    return {
        "status": "completed",
        "results": results,
        "summary": {
            "total_caches": len(results),
            "successful": successful,
            "failed": failed
        }
    }


@router.get("/health/system")
def get_system_health(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Get extended system health metrics (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s requested system health metrics", current_user.id)
    return get_extended_system_info()


@router.get("/health/queues")
def get_queue_health(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Get Celery queue statistics (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s requested queue health metrics", current_user.id)
    return get_celery_queue_stats()


@router.get("/health/connections")
def get_connection_health(current_user=Depends(get_current_user)) -> Dict[str, Any]:
    """Get database/Redis connection status (admin only)."""
    _require_admin(current_user)
    
    logger.info("Admin user %s requested connection health metrics", current_user.id)
    return get_connection_status()
