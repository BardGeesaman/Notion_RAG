"""
Email cleanup utilities.

This module handles cleanup operations for email/note ingestion:
- Finding and deleting orphaned chunks
- Deleting emails and their associated chunks
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

REQUESTS_PER_SECOND = 2.0

__all__ = [
    "cleanup_orphaned_chunks",
    "delete_email_and_chunks",
]


def _notion_base_url() -> str:
    """
    DEPRECATED: Notion support has been removed.
    Returns empty string for backward compatibility.
    """
    logger.debug("[EMAIL-CLEANUP] _notion_base_url() deprecated - Notion support removed")
    return ""


def _rag_db_id() -> str:
    """
    DEPRECATED: Notion support has been removed.
    Raises RuntimeError since Notion is no longer available.
    """
    logger.warning("[EMAIL-CLEANUP] _rag_db_id() deprecated - Notion support removed")
    raise RuntimeError("Notion support has been removed. Use Postgres-based cleanup instead.")


def cleanup_orphaned_chunks() -> None:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously found and deleted chunks whose parent email/note no longer exists in Notion.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Implement Postgres-based orphaned chunk cleanup if needed.
    """
    logger.warning(
        "[INGEST][EMAIL] cleanup_orphaned_chunks() deprecated - Notion support removed. "
        "Postgres is now the source of truth."
    )
    return


def delete_email_and_chunks(email_page_id: str) -> None:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously deleted an email and all its associated RAG chunks from Notion and Pinecone.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Implement Postgres-based email/chunk deletion if needed.
    
    Args:
        email_page_id: Notion page ID (ignored)
    """
    logger.warning(
        "[INGEST][EMAIL] delete_email_and_chunks() deprecated - Notion support removed. "
        "Postgres is now the source of truth."
    )
    return

