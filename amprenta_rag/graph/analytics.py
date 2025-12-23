"""Graph analytics utilities (NetworkX).

Computes lightweight metrics for visualization and ranking:
- Degree centrality
- Community detection (Louvain) + modularity
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Tuple

import networkx as nx
from networkx.algorithms import community


def _node_id(n: Any) -> Optional[str]:
    """Normalize node identifier to a stable string."""
    if n is None:
        return None
    if isinstance(n, str):
        return n
    if isinstance(n, dict):
        # Accept either explicit "id" or entity tuple fields
        if n.get("id"):
            return str(n["id"])
        et = n.get("entity_type")
        eid = n.get("entity_id")
        if et and eid:
            return f"{et}:{eid}"
        data = n.get("data") or {}
        if isinstance(data, dict) and data.get("id"):
            return str(data["id"])
    return None


def _edge_endpoints(e: Any) -> Optional[Tuple[str, str, float]]:
    """Normalize edge to (u, v, weight)."""
    if e is None:
        return None
    if isinstance(e, dict):
        # Prefer explicit endpoints
        if e.get("source") and e.get("target"):
            w = e.get("weight") or e.get("confidence") or 1.0
            try:
                return str(e["source"]), str(e["target"]), float(w)
            except Exception:
                return str(e["source"]), str(e["target"]), 1.0

        # GraphEdge-like dict
        st = e.get("source_entity_type")
        sid = e.get("source_entity_id")
        tt = e.get("target_entity_type")
        tid = e.get("target_entity_id")
        if st and sid and tt and tid:
            u = f"{st}:{sid}"
            v = f"{tt}:{tid}"
            w = e.get("confidence") or 1.0
            try:
                return u, v, float(w)
            except Exception:
                return u, v, 1.0

        # Cytoscape edge
        data = e.get("data") or {}
        if isinstance(data, dict) and data.get("source") and data.get("target"):
            w = data.get("confidence") or data.get("weight") or 1.0
            try:
                return str(data["source"]), str(data["target"]), float(w)
            except Exception:
                return str(data["source"]), str(data["target"]), 1.0
    return None


def compute_graph_analytics(
    nodes: Sequence[Any],
    edges: Sequence[Any],
    metrics: Sequence[str],
) -> Dict[str, Any]:
    """Compute degree centrality and/or community detection.

    Args:
        nodes: list of nodes; accepted formats:
            - {"id": "..."} or {"entity_type": "...", "entity_id": "..."}
            - Cytoscape nodes: {"data": {"id": "..."}}
        edges: list of edges; accepted formats:
            - {"source": "...", "target": "...", "confidence": 0.8}
            - GraphEdge dict with source_entity_type/source_entity_id/...
            - Cytoscape edges: {"data": {"source": "...", "target": "...", "confidence": ...}}
        metrics: list containing "degree_centrality" and/or "communities"

    Returns:
        dict with optional keys:
            - degree_centrality: {node_id: score}
            - communities: {node_id: community_id}
            - modularity: float
    """
    want_degree = "degree_centrality" in set(metrics or [])
    want_comms = "communities" in set(metrics or [])

    G = nx.Graph()

    for n in nodes or []:
        nid = _node_id(n)
        if nid:
            G.add_node(nid)

    for e in edges or []:
        parsed = _edge_endpoints(e)
        if not parsed:
            continue
        u, v, w = parsed
        if u not in G:
            G.add_node(u)
        if v not in G:
            G.add_node(v)
        G.add_edge(u, v, weight=float(w))

    out: Dict[str, Any] = {}

    if want_degree:
        out["degree_centrality"] = {k: float(v) for k, v in nx.degree_centrality(G).items()}

    if want_comms:
        if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
            out["communities"] = {}
            out["modularity"] = 0.0
        else:
            comms = community.louvain_communities(G, weight="weight", seed=42)
            node_to_comm: Dict[str, int] = {}
            for i, members in enumerate(comms):
                for nid in members:
                    node_to_comm[str(nid)] = int(i)
            out["communities"] = node_to_comm
            out["modularity"] = float(community.modularity(G, comms, weight="weight"))

    return out


__all__ = ["compute_graph_analytics"]


