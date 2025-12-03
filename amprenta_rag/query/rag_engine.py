# amprenta_rag/query/rag_engine.py

"""
High-level RAG orchestration.

This module handles:
- Match summarization and filtering
- Notion chunk retrieval
- Answer synthesis
- Complete RAG query orchestration
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from amprenta_rag.clients.notion_client import get_page_text
from amprenta_rag.clients.openai_client import (get_default_models,
                                                get_openai_client)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.pinecone_query import build_meta_filter, query_pinecone

logger = get_logger(__name__)


@dataclass
class MatchSummary:
    """Summary of a single match from Pinecone."""

    id: str
    score: float
    source: str
    title: str
    snippet: str
    tags: List[str]
    metadata: Dict[str, Any]


@dataclass
class RAGQueryResult:
    """Complete result from a RAG query."""

    query: str
    matches: List[MatchSummary]
    filtered_matches: List[MatchSummary]
    context_chunks: List[str]
    answer: str


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


def synthesize_answer(user_query: str, chunks: List[str]) -> str:
    """
    Use OpenAI to synthesize an answer given the query + context chunks.

    Args:
        user_query: User's query text
        chunks: List of context chunk texts from Notion

    Returns:
        Synthesized answer string from OpenAI
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()

    context = "\n\n---\n\n".join(chunks[:8])  # cap context for now

    system_prompt = (
        "You are an assistant helping interpret Amprenta's internal RAG context.\n"
        "Use only the provided context chunks to answer the user's question.\n"
        "If the context is insufficient, say so explicitly."
    )

    user_content = (
        f"User question:\n{user_query}\n\n" f"Relevant context chunks:\n{context}"
    )

    try:
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            temperature=0.2,
        )
    except Exception as e:
        logger.error("[RAG] OpenAI API error synthesizing answer: %r", e)
        raise
    return resp.choices[0].message.content.strip()  # type: ignore[union-attr]


def query_rag(
    user_query: str,
    *,
    top_k: int = 10,
    source_types: Optional[List[str]] = None,
    tag: Optional[str] = None,
    generate_answer: bool = True,
    disease: Optional[str] = None,
    target: Optional[str] = None,
    lipid: Optional[str] = None,
    signature: Optional[str] = None,
) -> RAGQueryResult:
    """
    High-level API: run a complete RAG query and get structured results.

    This function orchestrates the full RAG pipeline:
    1. Builds metadata filter from query constraints
    2. Queries Pinecone for matches (with source type filtering)
    3. Summarizes and filters matches
    4. Retrieves chunk text from Notion
    5. Synthesizes answer using OpenAI (if enabled)

    Args:
        user_query: User's query text
        top_k: Number of results to retrieve from Pinecone
        source_types: Optional list of source types to filter by (e.g., ["Literature", "Experiment"])
        tag: Optional tag substring filter
        generate_answer: Whether to synthesize an answer (default: True)
        disease: Optional disease filter
        target: Optional molecular target filter
        lipid: Optional lipid filter (canonical ID or raw label)
        signature: Optional lipid signature filter

    Returns:
        RAGQueryResult with matches, filtered matches, context chunks, and answer
    """
    logger.info("[RAG] Querying Pinecone (top_k=%d) for: %s", top_k, user_query)
    meta_filter = build_meta_filter(
        disease=disease,
        target=target,
        lipid=lipid,
        signature=signature,
    )

    raw_matches = query_pinecone(
        user_query,
        top_k=top_k,
        meta_filter=meta_filter,
        source_types=source_types,
    )
    if not raw_matches:
        return RAGQueryResult(
            query=user_query,
            matches=[],
            filtered_matches=[],
            context_chunks=[],
            answer="No matches returned from Pinecone.",
        )

    matches = [summarize_match(m) for m in raw_matches]
    filtered = filter_matches(matches, tag=tag)

    if not filtered:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=[],
            context_chunks=[],
            answer="No matches left after filtering (tag).",
        )

    chunks = collect_chunks(filtered)
    if not chunks:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=filtered,
            context_chunks=[],
            answer="No Notion chunks could be retrieved for the filtered matches.",
        )

    if not generate_answer:
        return RAGQueryResult(
            query=user_query,
            matches=matches,
            filtered_matches=filtered,
            context_chunks=chunks,
            answer="Answer generation disabled (generate_answer=False).",
        )

    answer = synthesize_answer(user_query, chunks)

    # Add provenance summary
    source_counts: Dict[str, int] = {}
    for m in filtered:
        src = m.metadata.get("source_type") or m.metadata.get("source") or "Unknown"
        source_counts[src] = source_counts.get(src, 0) + 1

    provenance_lines = [
        "\nSources used:",
        *[f"- {src}: {cnt} chunk(s)" for src, cnt in sorted(source_counts.items())],
    ]

    answer_with_provenance = answer + "\n\n" + "\n".join(provenance_lines)

    return RAGQueryResult(
        query=user_query,
        matches=matches,
        filtered_matches=filtered,
        context_chunks=chunks,
        answer=answer_with_provenance,
    )
