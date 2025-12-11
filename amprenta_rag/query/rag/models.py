"""
Data models for RAG query results.

This module defines dataclasses for structured RAG query results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List


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
    citations: List[Any] = field(default_factory=list)

