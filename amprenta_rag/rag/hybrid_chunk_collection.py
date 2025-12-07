"""
Hybrid chunk collection for Postgres + Notion RAG.

Collects chunks from both Postgres (structured data) and Notion (narrative)
for comprehensive RAG context.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Optional, Any

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.rag.postgres_resolver import (
    resolve_postgres_id_from_metadata,
    get_entity_type_from_metadata,
    fetch_postgres_context,
    get_notion_id_from_postgres,
)

# Avoid circular import: only import MatchSummary for type checking
if TYPE_CHECKING:
    from amprenta_rag.query.rag.models import MatchSummary

logger = get_logger(__name__)


def collect_hybrid_chunks(
    matches: List["MatchSummary"],
    prefer_postgres: bool = True,
) -> List[str]:
    """
    Collect chunks from matches, preferring Postgres data when available.
    
    This function:
    1. Checks for Postgres IDs in metadata
    2. Fetches structured data from Postgres
    3. Falls back to Notion if Postgres data unavailable
    4. Optionally combines both for richer context
    
    Args:
        matches: List of match summaries from Pinecone
        prefer_postgres: If True, use Postgres data when available
        
    Returns:
        List of chunk text strings
    """
    chunks = []
    
    for match in matches:
        metadata = match.metadata
        
        # Try to resolve Postgres ID
        postgres_id = resolve_postgres_id_from_metadata(metadata)
        entity_type = get_entity_type_from_metadata(metadata)
        
        if postgres_id and entity_type and prefer_postgres:
            # Fetch from Postgres
            postgres_text = fetch_postgres_context(postgres_id, entity_type)
            if postgres_text:
                chunks.append(postgres_text)
                logger.debug(
                    "[RAG][HYBRID] Collected chunk from Postgres: %s %s",
                    entity_type,
                    postgres_id,
                )
                continue
        
        # Fall back to Notion (existing behavior)
        # Check for notion_page_id in metadata
        notion_id = metadata.get("notion_page_id") or metadata.get("dataset_page_id")
        if notion_id:
            try:
                from amprenta_rag.ingestion.notion_pages import extract_page_content
                notion_content = extract_page_content(notion_id)
                if notion_content:
                    chunks.append(notion_content)
                    logger.debug(
                        "[RAG][HYBRID] Collected chunk from Notion: %s",
                        notion_id,
                    )
                    continue
            except Exception as e:
                logger.debug(
                    "[RAG][HYBRID] Could not fetch Notion content for %s: %r",
                    notion_id,
                    e,
                )
        
        # Last resort: use snippet from metadata
        snippet = metadata.get("snippet") or match.snippet
        if snippet:
            chunks.append(snippet)
            logger.debug("[RAG][HYBRID] Used snippet from metadata")
    
    return chunks


def collect_enhanced_chunks(
    matches: List["MatchSummary"],
    include_notion_narrative: bool = True,
) -> List[str]:
    """
    Collect enhanced chunks combining Postgres structured data with Notion narrative.
    
    This provides the richest context by combining:
    - Structured data from Postgres (relationships, metadata)
    - Narrative content from Notion (results, conclusions, descriptions)
    
    Args:
        matches: List of match summaries
        include_notion_narrative: Include Notion narrative if available
        
    Returns:
        List of enhanced chunk text strings
    """
    chunks = []
    
    for match in matches:
        metadata = match.metadata
        
        # Get Postgres data
        postgres_id = resolve_postgres_id_from_metadata(metadata)
        entity_type = get_entity_type_from_metadata(metadata)
        
        if postgres_id and entity_type:
            # Build text from Postgres (includes structured data)
            postgres_text = fetch_postgres_context(
                postgres_id,
                entity_type,
                include_notion_narrative=include_notion_narrative,
            )
            
            if postgres_text:
                chunks.append(postgres_text)
                continue
        
        # Fallback to Notion or snippet
        notion_id = metadata.get("notion_page_id") or metadata.get("dataset_page_id")
        if notion_id:
            try:
                from amprenta_rag.ingestion.notion_pages import extract_page_content
                notion_content = extract_page_content(notion_id)
                if notion_content:
                    chunks.append(notion_content)
                    continue
            except Exception:
                pass
        
        # Use snippet as last resort
        snippet = metadata.get("snippet") or match.snippet
        if snippet:
            chunks.append(snippet)
    
    return chunks

