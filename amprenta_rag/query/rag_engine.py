"""
Compatibility wrapper for the refactored RAG engine.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from rag_engine.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.query.rag import (
    RAGQueryResult,
    MatchSummary,
    collect_chunks,
    filter_matches,
    query_rag,
    signature_similarity_query,
    summarize_match,
    synthesize_answer,
)

__all__ = [
    "MatchSummary",
    "RAGQueryResult",
    "summarize_match",
    "filter_matches",
    "collect_chunks",
    "synthesize_answer",
    "query_rag",
    "signature_similarity_query",
]
