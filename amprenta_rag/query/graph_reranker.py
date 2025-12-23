"""GraphRAG reranker (graph-based boosting for RAG matches).

This module applies a lightweight score boost to RAG matches when the query
entities and chunk entities are connected in the Evidence Graph.
"""

from __future__ import annotations

import re
from dataclasses import replace
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

from amprenta_rag.query.rag.models import MatchSummary


Entity = Tuple[str, UUID]  # (entity_type, entity_id)


UUID_RE = re.compile(r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b", re.I)
TYPED_ID_RE = re.compile(
    r"\b(?P<typ>compound|signature|feature|pathway|dataset|program|experiment)\s*:\s*(?P<id>[0-9a-f-]{36})\b",
    re.I,
)


def _normalize_entity_type(t: str) -> str:
    s = (t or "").strip().lower()
    # Normalize common source_type variants from metadata
    mapping = {
        "compounds": "compound",
        "compound": "compound",
        "signature": "signature",
        "signatures": "signature",
        "feature": "feature",
        "features": "feature",
        "pathway": "pathway",
        "pathways": "pathway",
        "dataset": "dataset",
        "datasets": "dataset",
        "program": "program",
        "programs": "program",
        "experiment": "experiment",
        "experiments": "experiment",
    }
    return mapping.get(s, s)


def _extract_query_entities(query: str) -> List[Entity]:
    """Extract entities from a user query.

    Preference order:
    1) Use notebook entity_extractor (LLM + DB) if available.
    2) Fallback: parse explicit typed UUIDs like "compound:<uuid>".
    """
    query = query or ""

    # 1) Try the notebook entity extractor (best-effort; may require LLM/DB).
    try:
        from amprenta_rag.notebook.entity_extractor import extract_entities

        out = extract_entities(query)  # type: ignore[call-arg]
        # Expected shape: dict with keys like dataset_ids, experiment_ids, compound_ids, campaign_ids
        entities: List[Entity] = []
        mapping = {
            "dataset_ids": "dataset",
            "experiment_ids": "experiment",
            "compound_ids": "compound",
            "campaign_ids": "campaign",
            "signature_ids": "signature",
            "feature_ids": "feature",
            "program_ids": "program",
        }
        if isinstance(out, dict):
            for k, et in mapping.items():
                ids = out.get(k) or []
                if isinstance(ids, list):
                    for s in ids:
                        try:
                            entities.append((_normalize_entity_type(et), UUID(str(s))))
                        except Exception:
                            continue
        if entities:
            return entities
    except Exception:
        pass

    # 2) Regex fallback: typed UUIDs
    ents: List[Entity] = []
    for m in TYPED_ID_RE.finditer(query):
        try:
            ents.append((_normalize_entity_type(m.group("typ")), UUID(m.group("id"))))
        except Exception:
            continue
    return ents


def _extract_chunk_entity(meta: Dict[str, Any]) -> Optional[Entity]:
    """Extract a single 'entity' identifier from chunk metadata."""
    if not isinstance(meta, dict):
        return None

    # Prefer RAGChunk fields
    st = meta.get("source_type") or meta.get("sourceType") or meta.get("entity_type") or meta.get("entityType")
    sid = meta.get("source_id") or meta.get("sourceId") or meta.get("entity_id") or meta.get("entityId")

    if st and sid:
        try:
            return (_normalize_entity_type(str(st)), UUID(str(sid)))
        except Exception:
            return None
    return None


@lru_cache(maxsize=20000)
def _cached_path_boost(
    src_type: str,
    src_id: str,
    tgt_type: str,
    tgt_id: str,
    max_hops: int,
    rels: Tuple[str, ...],
) -> float:
    """Compute a graph boost between two entities using shortest_path()."""
    try:
        from amprenta_rag.graph.traversal import shortest_path
    except Exception:
        return 0.0

    try:
        out = shortest_path(
            src_type,
            UUID(src_id),
            tgt_type,
            UUID(tgt_id),
            relationships=list(rels) if rels else None,
            timeout_s=5.0,
            max_nodes=500,
        )
    except Exception:
        return 0.0

    if not out:
        return 0.0

    nodes = out.get("nodes") or []
    edges = out.get("edges") or []
    try:
        path_len = max(0, len(nodes) - 1)
    except Exception:
        return 0.0

    if path_len == 0 or path_len > int(max_hops):
        return 0.0

    confs: List[float] = []
    for e in edges:
        c = (e or {}).get("confidence")
        if c is None:
            continue
        try:
            confs.append(float(c))
        except Exception:
            continue
    path_conf = min(confs) if confs else 0.5

    # Spec: graph_boost = confidence * 0.5^(path_length - 1)
    decay = 0.5 ** max(0, path_len - 1)
    return max(0.0, float(path_conf) * float(decay))


def _clear_path_cache() -> None:
    _cached_path_boost.cache_clear()


def boost_by_graph(
    query: str,
    matches: List[MatchSummary],
    alpha: float = 0.3,
    max_hops: int = 2,
) -> List[MatchSummary]:
    """Apply graph-based boosting to RAG matches.

    Algorithm:
    1) Extract entities from query
    2) Extract chunk entity from match.metadata
    3) For each (query_entity, chunk_entity) compute:
       graph_boost = confidence * 0.5^(path_length - 1), only if shortest path <= max_hops
    4) Blend:
       final_score = score * (1 + alpha * graph_boost)
    5) Return matches sorted by final_score desc
    """
    if not matches:
        return matches

    try:
        alpha_f = float(alpha)
    except Exception:
        alpha_f = 0.3

    q_ents = _extract_query_entities(query)
    if not q_ents:
        return matches

    boosted: List[Tuple[int, MatchSummary]] = []
    rels: Tuple[str, ...] = tuple()

    for idx, m in enumerate(matches):
        meta = m.metadata or {}
        c_ent = _extract_chunk_entity(meta)
        if c_ent is None:
            boosted.append((idx, m))
            continue

        best_boost = 0.0
        for (qt, qid) in q_ents:
            b = _cached_path_boost(qt, str(qid), c_ent[0], str(c_ent[1]), int(max_hops), rels)
            if b > best_boost:
                best_boost = b

        if best_boost <= 0.0:
            boosted.append((idx, m))
            continue

        new_score = float(m.score) * (1.0 + alpha_f * best_boost)
        new_meta = dict(meta)
        new_meta["graph_boost"] = best_boost
        new_meta["graph_alpha"] = alpha_f
        boosted.append((idx, replace(m, score=new_score, metadata=new_meta)))

    boosted.sort(key=lambda it: (-float(it[1].score), it[0]))
    return [m for _, m in boosted]


__all__ = ["boost_by_graph", "_clear_path_cache"]


