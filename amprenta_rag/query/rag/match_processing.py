"""
Match processing and summarization for RAG queries.

This module provides functions for processing and filtering Pinecone matches.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.query.rag.models import MatchSummary


def _get_source(meta: Dict[str, Any]) -> str:
    """Extract source from metadata."""
    return meta.get("source") or meta.get("source_type") or "Unknown"


def _format_snippet(meta: Dict[str, Any]) -> str:
    """Format snippet from metadata."""
    snippet = (meta.get("snippet") or "").replace("\n", " ")
    return snippet[:200]


def _get_tags(meta: Dict[str, Any]) -> List[str]:
    """Extract tags from metadata."""
    tags = meta.get("tags") or meta.get("zotero_tags") or []
    if isinstance(tags, str):
        tags = [tags]
    return list(tags)


def summarize_match(raw_match: Any) -> MatchSummary:
    """
    Convert a raw Pinecone match into a MatchSummary.

    Args:
        raw_match: Raw match object from Pinecone (with metadata)

    Returns:
        MatchSummary object with extracted fields
    """
    meta: Dict[str, Any] = raw_match.metadata or raw_match.get("metadata", {})  # type: ignore
    source = _get_source(meta)
    title = meta.get("title") or "(untitled)"
    snippet = _format_snippet(meta)
    tags = _get_tags(meta)

    return MatchSummary(
        id=getattr(raw_match, "id", None) or raw_match.get("id"),
        score=getattr(raw_match, "score", None) or raw_match.get("score", 0.0),
        source=source,
        title=title,
        snippet=snippet,
        tags=tags,
        metadata=meta,
    )


def filter_matches(
    matches: List[MatchSummary],
    tag: Optional[str] = None,
) -> List[MatchSummary]:
    """
    Apply optional tag filter to matches.

    Note: Source type filtering is now handled at the Pinecone query level.

    Args:
        matches: List of MatchSummary objects to filter
        tag: Optional tag substring to filter by

    Returns:
        Filtered list of MatchSummary objects
    """

    def _has_tag(m: MatchSummary) -> bool:
        if not tag:
            return True
        return any(tag.lower() in t.lower() for t in m.tags)

    return [m for m in matches if _has_tag(m)]

