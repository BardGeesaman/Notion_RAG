"""System health and monitoring utilities."""
from __future__ import annotations

import platform
import sys
import time
from typing import Any, Dict

from sqlalchemy.orm import Session

from amprenta_rag.database.models import User, Experiment, Compound, Signature, Dataset
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_db_stats(db: Session) -> Dict[str, int]:
    """
    Get database statistics (row counts for key tables).

    Args:
        db: Database session

    Returns:
        Dict with counts for users, experiments, compounds, signatures, datasets
    """
    stats = {
        "users": db.query(User).filter(User.is_active.is_(True)).count(),
        "experiments": db.query(Experiment).count(),
        "compounds": db.query(Compound).count(),
        "signatures": db.query(Signature).count(),
        "datasets": db.query(Dataset).count(),
    }

    return stats


def check_api_connectivity() -> Dict[str, bool]:
    """
    Check connectivity to external APIs (OpenAI, Pinecone).

    Returns:
        Dict with connectivity status: {openai: bool, pinecone: bool}
    """
    result = {"openai": False, "pinecone": False}

    # Check OpenAI
    try:
        from amprenta_rag.clients.openai_client import get_openai_client
        openai_client: Any = get_openai_client()
        # Try a simple API call (list models is lightweight)
        openai_client.models.list()
        result["openai"] = True
    except Exception as e:
        logger.debug("[HEALTH] OpenAI connectivity check failed: %r", e)

    # Pinecone backend deprecated - using pgvector instead
    result["pinecone"] = False  # No longer supported

    return result


def get_system_info() -> Dict[str, Any]:
    """
    Get system information (Python version, platform, memory).

    Returns:
        Dict with system information
    """
    info: Dict[str, Any] = {
        "python_version": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "architecture": platform.machine(),
    }

    # Try to get memory info if psutil is available
    try:
        import psutil
        memory = psutil.virtual_memory()
        info["memory"] = {
            "total_gb": round(memory.total / (1024**3), 2),
            "available_gb": round(memory.available / (1024**3), 2),
            "used_gb": round(memory.used / (1024**3), 2),
            "percent": memory.percent,
        }
    except ImportError:
        info["memory"] = {"available": False, "note": "psutil not installed"}

    return info


def get_extended_system_info() -> Dict[str, Any]:
    """Get extended system metrics (CPU, Memory, Disk)."""
    try:
        import psutil
        
        # Get CPU information
        cpu_info = {
            "cpu_percent": psutil.cpu_percent(interval=0.1),
            "cpu_count": psutil.cpu_count(),
            "cpu_count_logical": psutil.cpu_count(logical=True),
        }
        
        # Get memory information
        memory = psutil.virtual_memory()
        memory_info = {
            "total_gb": round(memory.total / (1024**3), 2),
            "available_gb": round(memory.available / (1024**3), 2),
            "used_gb": round(memory.used / (1024**3), 2),
            "percent": memory.percent,
        }
        
        # Get disk information
        try:
            disk = psutil.disk_usage('/')
            disk_info = {
                "total_gb": round(disk.total / (1024**3), 2),
                "free_gb": round(disk.free / (1024**3), 2),
                "used_gb": round(disk.used / (1024**3), 2),
                "percent": round((disk.used / disk.total) * 100, 2),
            }
        except Exception as e:
            logger.warning("Failed to get disk usage: %s", e)
            disk_info = {"error": str(e), "status": "unavailable"}
        
        return {
            "cpu": cpu_info,
            "memory": memory_info,
            "disk": disk_info,
            "timestamp": time.time(),
        }
        
    except ImportError:
        logger.warning("psutil not available for extended system monitoring")
        return {
            "error": "psutil not installed",
            "status": "unavailable",
            "cpu": {"status": "unavailable"},
            "memory": {"status": "unavailable"},
            "disk": {"status": "unavailable"},
        }
    except Exception as e:
        logger.error("Failed to get extended system info: %s", e)
        return {
            "error": str(e),
            "status": "error",
        }


def get_celery_queue_stats() -> Dict[str, Any]:
    """Get Celery queue statistics from Redis."""
    try:
        from amprenta_rag.jobs.celery_app import celery_app
        
        # Get inspection interface
        inspect = celery_app.control.inspect()
        
        # Get various queue statistics
        active_tasks = inspect.active() or {}
        reserved_tasks = inspect.reserved() or {}
        scheduled_tasks = inspect.scheduled() or {}
        
        # Count total tasks across all workers
        total_active = sum(len(worker_tasks) for worker_tasks in active_tasks.values())
        total_reserved = sum(len(worker_tasks) for worker_tasks in reserved_tasks.values())
        total_scheduled = sum(len(worker_tasks) for worker_tasks in scheduled_tasks.values())
        
        return {
            "active": total_active,
            "reserved": total_reserved,
            "scheduled": total_scheduled,
            "workers": list(active_tasks.keys()) if active_tasks else [],
            "worker_count": len(active_tasks) if active_tasks else 0,
            "status": "available",
        }
        
    except ImportError as e:
        logger.warning("Celery not available: %s", e)
        return {
            "error": f"Celery not available: {str(e)}",
            "status": "unavailable",
        }
    except Exception as e:
        logger.error("Failed to get Celery queue stats: %s", e)
        return {
            "error": str(e),
            "status": "error",
        }


def get_connection_status() -> Dict[str, Any]:
    """Get database and Redis connection status."""
    status = {}
    
    # Check PostgreSQL connection
    try:
        from amprenta_rag.database.session import db_session
        from sqlalchemy import text
        with db_session() as db:
            # Simple query to test connection
            db.execute(text("SELECT 1"))
            status["postgresql"] = {
                "status": "connected",
                "type": "PostgreSQL",
            }
    except Exception as e:
        logger.error("PostgreSQL connection check failed: %s", e)
        status["postgresql"] = {
            "status": "error",
            "error": str(e),
            "type": "PostgreSQL",
        }
    
    # Check Redis connection (used by Celery)
    try:
        from amprenta_rag.jobs.celery_app import celery_app
        
        # Try to get broker connection info
        broker_url = celery_app.conf.broker_url
        if broker_url and "redis" in broker_url.lower():
            # Try to connect to Redis
            import redis
            # Parse Redis URL to get connection details
            if broker_url.startswith('redis://'):
                redis_client = redis.from_url(broker_url)
                redis_client.ping()
                status["redis"] = {
                    "status": "connected",
                    "type": "Redis",
                    "broker_url": broker_url.split('@')[-1] if '@' in broker_url else broker_url,  # Hide credentials
                }
            else:
                status["redis"] = {
                    "status": "unavailable",
                    "type": "Redis",
                    "note": "Non-Redis broker configured",
                }
        else:
            status["redis"] = {
                "status": "unavailable",
                "type": "Redis",
                "note": "Redis not configured as Celery broker",
            }
            
    except ImportError as e:
        logger.warning("Redis client not available: %s", e)
        status["redis"] = {
            "status": "unavailable",
            "error": f"Redis client not available: {str(e)}",
            "type": "Redis",
        }
    except Exception as e:
        logger.error("Redis connection check failed: %s", e)
        status["redis"] = {
            "status": "error",
            "error": str(e),
            "type": "Redis",
        }
    
    return status
