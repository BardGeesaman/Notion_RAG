"""BM25/Full-text search using PostgreSQL tsvector."""
from __future__ import annotations

from typing import List, Dict, Any, Optional
from sqlalchemy import func, text

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import RAGChunk
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def bm25_search(
    query: str,
    limit: int = 20,
    source_type: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Perform full-text search on RAG chunks using PostgreSQL.

    Args:
        query: Search query text
        limit: Max results to return
        source_type: Optional filter by source type

    Returns:
        List of dicts with chunk_id, content, score, metadata
    """
    if not query.strip():
        return []

    db_gen = get_db()
    db = next(db_gen)
    try:
        ts_query = func.plainto_tsquery('english', query)
        base = db.query(
            RAGChunk.id,
            RAGChunk.chunk_text,
            RAGChunk.source_type,
            RAGChunk.chunk_metadata,
            func.ts_rank(RAGChunk.search_vector, ts_query).label('rank')
        ).filter(
            RAGChunk.search_vector.op('@@')(ts_query)
        )
        if source_type:
            base = base.filter(RAGChunk.source_type == source_type)
        results = base.order_by(text('rank DESC')).limit(limit).all()
        return [
            {
                "chunk_id": str(r.id),
                "content": r.chunk_text,
                "source_type": r.source_type,
                "metadata": r.chunk_metadata or {},
                "bm25_score": float(r.rank) if r.rank else 0.0,
            }
            for r in results
        ]
    except Exception as e:
        logger.warning("[BM25] Search failed: %s", e)
        return []
    finally:
        db_gen.close()


def reciprocal_rank_fusion(
    bm25_results: List[Dict[str, Any]],
    vector_results: List[Dict[str, Any]],
    k: int = 60,
    alpha: float = 0.5,
) -> List[Dict[str, Any]]:
    """
    Combine BM25 and vector results using Reciprocal Rank Fusion.
    """
    scores = {}
    result_map = {}

    for rank, r in enumerate(bm25_results):
        chunk_id = r.get("chunk_id") or r.get("id")
        rrf_score = (1 - alpha) * (1 / (k + rank + 1))
        scores[chunk_id] = scores.get(chunk_id, 0) + rrf_score
        result_map[chunk_id] = r

    for rank, r in enumerate(vector_results):
        chunk_id = r.get("chunk_id") or r.get("id")
        rrf_score = alpha * (1 / (k + rank + 1))
        scores[chunk_id] = scores.get(chunk_id, 0) + rrf_score
        if chunk_id not in result_map:
            result_map[chunk_id] = r

    sorted_ids = sorted(scores.keys(), key=lambda x: scores[x], reverse=True)

    return [
        {**result_map[cid], "rrf_score": scores[cid]}
        for cid in sorted_ids
    ]
