"""Cytoscape JSON helpers for graph visualization."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Tuple
from uuid import UUID


def _node_key(entity_type: str, entity_id: UUID) -> str:
    return f"{entity_type}:{entity_id}"


def to_cytoscape_json(
    nodes: Iterable[Tuple[str, UUID]],
    edges: Iterable[Dict[str, Any]],
    node_labels: Mapping[Tuple[str, UUID], str] | None = None,
) -> Dict[str, Any]:
    """Convert node/edge lists to Cytoscape-compatible JSON.

    Args:
        nodes: iterable of (entity_type, entity_id)
        edges: iterable of edge dicts containing:
            source_entity_type, source_entity_id, target_entity_type, target_entity_id,
            relationship_type, confidence, evidence_source, provenance
    """
    cy_nodes: List[Dict[str, Any]] = []
    for et, eid in nodes:
        fallback = f"{et}:{str(eid)[:8]}"
        label = (node_labels or {}).get((et, eid), fallback)
        cy_nodes.append(
            {
                "data": {
                    "id": _node_key(et, eid),
                    "entity_type": et,
                    "entity_id": str(eid),
                    "label": label,
                }
            }
        )

    cy_edges: List[Dict[str, Any]] = []
    for e in edges:
        s_type = e["source_entity_type"]
        s_id = UUID(str(e["source_entity_id"]))
        t_type = e["target_entity_type"]
        t_id = UUID(str(e["target_entity_id"]))
        rel = e.get("relationship_type") or ""
        cy_edges.append(
            {
                "data": {
                    "id": e.get("id") or f"{_node_key(s_type, s_id)}->{_node_key(t_type, t_id)}:{rel}",
                    "source": _node_key(s_type, s_id),
                    "target": _node_key(t_type, t_id),
                    "relationship_type": rel,
                    "confidence": e.get("confidence"),
                    "evidence_source": e.get("evidence_source"),
                    "provenance": e.get("provenance"),
                }
            }
        )

    return {"nodes": cy_nodes, "edges": cy_edges}


__all__ = ["to_cytoscape_json"]


