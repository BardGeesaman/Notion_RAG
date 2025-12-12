"""
Data models for RAG query results.

This module defines dataclasses for structured RAG query results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from amprenta_rag.query.evaluation import EvalResult


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
class Citation:
    """Citation metadata for RAG answers."""

    number: int
    chunk_id: str
    source_type: str | None = None
    title: str | None = None
    dataset_name: str | None = None
    experiment_name: str | None = None
    url: str | None = None


@dataclass
class RAGQueryResult:
    """Complete result from a RAG query."""

    query: str
    matches: List[MatchSummary]
    filtered_matches: List[MatchSummary]
    context_chunks: List[str]
    answer: str
    citations: List[Citation] = field(default_factory=list)
    groundedness_score: Optional[float] = None
    unsupported_claims: List[str] = field(default_factory=list)
    evaluation: Optional["EvalResult"] = None
    trust_summary: Optional[Dict[str, Any]] = None

