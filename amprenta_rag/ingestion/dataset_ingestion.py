# amprenta_rag/ingestion/dataset_ingestion.py
"""
Dataset ingestion module.

Handles ingestion of experimental datasets from Notion into Pinecone.
Supports mwTab data extraction, species extraction, signature matching,
and automatic metadata population.
"""

from __future__ import annotations

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def ingest_dataset(page_id: str, force: bool = False) -> None:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously ingested datasets from Notion.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Remove this function or refactor to use Postgres UUIDs.
    
    Args:
        page_id: Notion page ID (ignored)
        force: If True, re-ingest (ignored)
    """
    logger.warning(
        "[INGEST][DATASET] ingest_dataset() deprecated - Notion support removed. "
        "Use Postgres-based ingestion instead."
    )
    return
