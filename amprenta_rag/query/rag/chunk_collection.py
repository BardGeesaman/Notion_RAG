"""
Chunk collection from Notion for RAG queries.

This module provides functions for retrieving full chunk text from Notion.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.clients.notion_client import get_page_text
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.rag.models import MatchSummary

logger = get_logger(__name__)


def collect_chunks(matches: List[MatchSummary]) -> List[str]:
    """
    Retrieve full chunk text from Notion for each match.

    Uses the 'notion_chunk_page_id' metadata key to fetch chunk text.

    Args:
        matches: List of MatchSummary objects with notion_chunk_page_id in metadata

    Returns:
        List of chunk text strings from Notion
    """
    chunks: List[str] = []
    for m in matches:
        page_id = m.metadata.get("notion_chunk_page_id")
        if not page_id:
            continue
        try:
            text = get_page_text(page_id)
        except Exception as e:
            logger.warning("[RAG] Failed to fetch Notion chunk %s: %s", page_id, e)
            continue
        if text:
            chunks.append(text)
    return chunks

