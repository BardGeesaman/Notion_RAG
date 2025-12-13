"""
RAG query engine.

This package provides functions for querying the RAG system, including
match processing, chunk collection, answer synthesis, and signature similarity queries.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.query.rag.chunk_collection import collect_chunks
from amprenta_rag.query.rag.match_processing import filter_matches, summarize_match
from amprenta_rag.query.rag.models import MatchSummary, RAGQueryResult
from amprenta_rag.query.rag.query import query_rag
from amprenta_rag.query.rag.synthesis import synthesize_answer

__all__ = [
    "MatchSummary",
    "RAGQueryResult",
    "summarize_match",
    "filter_matches",
    "collect_chunks",
    "synthesize_answer",
    "query_rag",
]

