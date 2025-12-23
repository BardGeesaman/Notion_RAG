"""Graph traversal utilities for Evidence Graph.

Constraints:
- Max 500 nodes per response.
- Path finding timeout: 5s.
"""

from __future__ import annotations

import time
from collections import deque
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
from uuid import UUID

from amprenta_rag.graph.edge_builder import EdgeBuilder


Node = Tuple[str, UUID]
EdgeDict = Dict[str, Any]


def _edge_to_dict(n) -> EdgeDict:
    return {
        "source_entity_type": n.source_entity_type,
        "source_entity_id": n.source_entity_id,
        "target_entity_type": n.target_entity_type,
        "target_entity_id": n.target_entity_id,
        "relationship_type": n.relationship_type,
        "confidence": n.confidence,
        "evidence_source": n.evidence_source,
        "provenance": n.provenance,
    }


def k_hop_subgraph(
    entity_type: str,
    entity_id: UUID,
    depth: int = 2,
    relationships: Optional[Iterable[str]] = None,
    *,
    max_nodes: int = 500,
) -> Dict[str, Any]:
    """BFS expansion up to depth hops returning nodes+edges."""
    depth = max(0, int(depth))
    builder = EdgeBuilder()

    visited: Set[Node] = set()
    visited_edges: Set[Tuple[str, str, str, str, str]] = set()
    nodes: List[Node] = []
    edges: List[EdgeDict] = []

    start: Node = (entity_type, entity_id)
    frontier: List[Tuple[Node, int]] = [(start, 0)]

    while frontier:
        (n_type, n_id), dist = frontier.pop(0)
        node = (n_type, n_id)
        if node in visited:
            continue
        visited.add(node)
        nodes.append(node)
        if len(nodes) >= max_nodes:
            break

        if dist >= depth:
            continue

        neighs = builder.get_neighbors(
            entity_type=n_type,
            entity_id=n_id,
            direction="both",
            relationship_types=relationships,
            limit=2000,
        )

        for nb in neighs:
            ekey = (
                nb.source_entity_type,
                str(nb.source_entity_id),
                nb.target_entity_type,
                str(nb.target_entity_id),
                nb.relationship_type,
            )
            if ekey not in visited_edges:
                visited_edges.add(ekey)
                edges.append(_edge_to_dict(nb))

            # Add the opposite endpoint as next node
            if (nb.source_entity_type, nb.source_entity_id) == node:
                other = (nb.target_entity_type, nb.target_entity_id)
            else:
                other = (nb.source_entity_type, nb.source_entity_id)

            if other not in visited and len(nodes) + len(frontier) < max_nodes:
                frontier.append((other, dist + 1))

    return {"nodes": nodes, "edges": edges, "truncated": len(nodes) >= max_nodes}


def shortest_path(
    source_type: str,
    source_id: UUID,
    target_type: str,
    target_id: UUID,
    relationships: Optional[Iterable[str]] = None,
    *,
    timeout_s: float = 5.0,
    max_nodes: int = 500,
) -> Optional[Dict[str, Any]]:
    """Bidirectional BFS shortest path. Returns None if no path or timeout."""
    start_t = time.monotonic()
    builder = EdgeBuilder()

    src: Node = (source_type, source_id)
    tgt: Node = (target_type, target_id)
    if src == tgt:
        return {"nodes": [src], "edges": []}

    # Frontiers + parent maps
    q_f = deque([src])
    q_b = deque([tgt])
    parent_f: Dict[Node, Tuple[Optional[Node], Optional[EdgeDict]]] = {src: (None, None)}
    parent_b: Dict[Node, Tuple[Optional[Node], Optional[EdgeDict]]] = {tgt: (None, None)}

    def _expand_one(
        q: deque,
        parents_this: Dict[Node, Tuple[Optional[Node], Optional[EdgeDict]]],
        parents_other: Dict[Node, Tuple[Optional[Node], Optional[EdgeDict]]],
    ) -> Optional[Node]:
        if not q:
            return None
        cur = q.popleft()
        neighs = builder.get_neighbors(
            entity_type=cur[0],
            entity_id=cur[1],
            direction="both",
            relationship_types=relationships,
            limit=2000,
        )
        for nb in neighs:
            e = _edge_to_dict(nb)
            # Determine neighbor node
            if (nb.source_entity_type, nb.source_entity_id) == cur:
                nxt = (nb.target_entity_type, nb.target_entity_id)
            else:
                nxt = (nb.source_entity_type, nb.source_entity_id)

            if nxt in parents_this:
                continue
            parents_this[nxt] = (cur, e)
            if nxt in parents_other:
                return nxt
            if len(parents_this) + len(parents_other) >= max_nodes:
                return None
            q.append(nxt)
        return None

    meet: Optional[Node] = None
    while q_f and q_b:
        if time.monotonic() - start_t > timeout_s:
            return None

        # Expand smaller frontier
        if len(q_f) <= len(q_b):
            meet = _expand_one(q_f, parent_f, parent_b)
        else:
            meet = _expand_one(q_b, parent_b, parent_f)
        if meet is not None:
            break

    if meet is None:
        return None

    # Reconstruct path nodes + edges
    path_nodes_f: List[Node] = []
    path_edges_f: List[EdgeDict] = []
    cur = meet
    while cur is not None:
        prev, e = parent_f[cur]
        path_nodes_f.append(cur)
        if e is not None:
            path_edges_f.append(e)
        cur = prev
    path_nodes_f.reverse()
    path_edges_f.reverse()

    path_nodes_b: List[Node] = []
    path_edges_b: List[EdgeDict] = []
    cur = meet
    prev_b, _ = parent_b[cur]
    cur = prev_b
    while cur is not None:
        prev, e = parent_b[cur]
        path_nodes_b.append(cur)
        if e is not None:
            path_edges_b.append(e)
        cur = prev

    nodes = path_nodes_f + path_nodes_b
    edges = path_edges_f + path_edges_b
    return {"nodes": nodes, "edges": edges}


__all__ = ["k_hop_subgraph", "shortest_path"]


