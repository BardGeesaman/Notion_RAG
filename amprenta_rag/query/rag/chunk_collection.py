"""
Chunk collection from Notion for RAG queries.

This module provides functions for retrieving full chunk text from Notion.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.rag.models import MatchSummary

logger = get_logger(__name__)


def collect_chunks(matches: List[MatchSummary]) -> List[str]:
    """
    DEPRECATED: Notion support has been removed.
    
    This function previously retrieved full chunk text from Notion.
    Postgres is now the source of truth.
    
    TODO: Phase 3 - Implement Postgres-based chunk retrieval or remove this function.
    
    Args:
        matches: List of MatchSummary objects

    Returns:
        Empty list (Notion support removed)
    """
    logger.debug(
        "[RAG][CHUNK-COLLECTION] collect_chunks() deprecated - Notion support removed"
    )
    # Return empty chunks - use snippet from metadata if available
    chunks: List[str] = []
    for m in matches:
        snippet = m.metadata.get("snippet") or m.snippet
        if snippet:
            chunks.append(snippet)
    return chunks

