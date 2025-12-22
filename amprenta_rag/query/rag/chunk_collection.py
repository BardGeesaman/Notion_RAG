"""
Chunk collection for RAG queries.

Historically, `collect_chunks()` returned a list of context strings built from
the `snippet` fields in match results.

In the Postgres-backed system, we can additionally fetch chunk text from the
`RAGChunk` table when a match provides a `chunk_id`.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.rag.models import MatchSummary
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import RAGChunk

logger = get_logger(__name__)


def _get_match_chunk_id(match: MatchSummary) -> Optional[str]:
    chunk_id = getattr(match, "chunk_id", None)
    if chunk_id:
        return str(chunk_id)
    # Back-compat: sometimes the match id itself is the chunk id
    if getattr(match, "id", None):
        return str(match.id)
    return None


def collect_chunks(matches: List[MatchSummary]) -> List[str]:
    """
    Collect context chunks for the given matches.

    Priority:
    - Use `match.snippet` if present (fast, no DB required).
    - Otherwise, attempt to load the full chunk text from Postgres via `RAGChunk`.
    """
    # 1) Use snippets where available
    snippets: List[Optional[str]] = []
    missing_ids: List[str] = []
    for m in matches:
        snippet = getattr(m, "snippet", None)
        if isinstance(snippet, str) and snippet.strip():
            snippets.append(snippet)
            continue
        snippets.append(None)
        mid = _get_match_chunk_id(m)
        if mid:
            missing_ids.append(mid)

    id_to_text: Dict[str, str] = {}
    if missing_ids:
        try:
            with db_session() as db:
                rows = db.query(RAGChunk).filter(RAGChunk.id.in_(missing_ids)).all()
                for row in rows:
                    text = getattr(row, "chunk_text", None) or getattr(row, "chunk", None)
                    if text:
                        id_to_text[str(row.id)] = str(text)
        except Exception as e:
            logger.warning("[RAG] Failed to fetch chunks from Postgres: %r", e)

    # 2) Return in the original order, skipping missing content
    out: List[str] = []
    for idx, m in enumerate(matches):
        snippet_val = snippets[idx]
        if snippet_val:
            out.append(snippet_val)
            continue
        mid = _get_match_chunk_id(m)
        if mid and mid in id_to_text:
            out.append(id_to_text[mid])
    return out

