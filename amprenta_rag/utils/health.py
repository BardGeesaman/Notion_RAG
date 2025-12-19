"""System health and monitoring utilities."""
from __future__ import annotations

import platform
import sys
from typing import Dict, Any

from amprenta_rag.database.models import User, Experiment, Compound, Signature, Dataset
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_db_stats(db) -> Dict[str, int]:
    """
    Get database statistics (row counts for key tables).

    Args:
        db: Database session

    Returns:
        Dict with counts for users, experiments, compounds, signatures, datasets
    """
    stats = {
        "users": db.query(User).filter(User.is_active).count(),
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
        client = get_openai_client()
        # Try a simple API call (list models is lightweight)
        client.models.list()
        result["openai"] = True
    except Exception as e:
        logger.debug("[HEALTH] OpenAI connectivity check failed: %r", e)

    # Check Pinecone
    try:
        from amprenta_rag.clients.pinecone_client import get_pinecone_client
        client = get_pinecone_client()
        # Try listing indexes (lightweight operation)
        client.list_indexes()
        result["pinecone"] = True
    except Exception as e:
        logger.debug("[HEALTH] Pinecone connectivity check failed: %r", e)

    return result


def get_system_info() -> Dict[str, Any]:
    """
    Get system information (Python version, platform, memory).

    Returns:
        Dict with system information
    """
    info = {
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
