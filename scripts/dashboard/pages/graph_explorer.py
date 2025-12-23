"""Graph Explorer dashboard page (Evidence Graph)."""

from __future__ import annotations

import os
from typing import Any, Dict
from uuid import UUID

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")
TRAVERSE_ENDPOINT = f"{API_BASE}/api/graph/traverse"
PATH_ENDPOINT = f"{API_BASE}/api/graph/path"
ANALYTICS_ENDPOINT = f"{API_BASE}/api/graph/analytics"


def _to_agraph(
    cytoscape: Dict[str, Any],
    *,
    degree_centrality: Dict[str, float] | None = None,
    communities: Dict[str, int] | None = None,
    size_by_degree: bool = False,
    color_by_community: bool = False,
):
    from streamlit_agraph import Node, Edge  # type: ignore[import-not-found]

    COLORS = {
        "compound": "#2E86AB",
        "signature": "#A23B72",
        "feature": "#F18F01",
        "pathway": "#C73E1D",
        "dataset": "#3B1F2B",
    }

    # Simple palette for community recoloring (fallback if > len palette)
    COMM_COLORS = [
        "#2E86AB",
        "#A23B72",
        "#F18F01",
        "#C73E1D",
        "#3B1F2B",
        "#4ECDC4",
        "#C7F9CC",
        "#9B5DE5",
        "#F15BB5",
        "#00BBF9",
    ]

    nodes = []
    deg = degree_centrality or {}
    comms = communities or {}
    max_deg = max(deg.values()) if deg else 0.0
    for n in cytoscape.get("nodes", []):
        d = (n or {}).get("data") or {}
        base_color = COLORS.get(d.get("entity_type"), "#888888")
        node_id = d.get("id")
        color = base_color
        if color_by_community and node_id and node_id in comms:
            cid = int(comms.get(node_id, 0))
            color = COMM_COLORS[cid % len(COMM_COLORS)]
        size = 20
        if size_by_degree and node_id and max_deg > 0:
            val = float(deg.get(node_id, 0.0))
            # Map [0, max] -> [10, 50]
            size = int(10 + (40 * (val / max_deg)))
        nodes.append(
            Node(
                id=d.get("id"),
                label=d.get("label"),
                size=size,
                color=color,
                title=f"{d.get('entity_type')} {d.get('entity_id')}",
            )
        )
    edges = []
    for e in cytoscape.get("edges", []):
        d = (e or {}).get("data") or {}
        edges.append(
            Edge(
                source=d.get("source"),
                target=d.get("target"),
                label=d.get("relationship_type"),
            )
        )
    return nodes, edges


def render_graph_explorer_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🕸️ Graph Explorer")
    st.caption("Explore evidence graph neighborhoods and shortest paths.")

    # Deep-link support: /?page=Graph%20Explorer&entity_type=...&entity_id=...
    qp = getattr(st, "query_params", {})
    default_type = None
    default_id = None
    try:
        if "entity_type" in qp:
            default_type = qp["entity_type"][0] if isinstance(qp["entity_type"], list) else str(qp["entity_type"])
        if "entity_id" in qp:
            raw = qp["entity_id"][0] if isinstance(qp["entity_id"], list) else str(qp["entity_id"])
            default_id = raw
    except Exception:
        pass

    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        entity_type = st.selectbox(
            "Entity type",
            ["compound", "signature", "feature", "pathway", "dataset", "program", "experiment"],
            index=0 if default_type is None else max(0, ["compound", "signature", "feature", "pathway", "dataset", "program", "experiment"].index(default_type))  # type: ignore[arg-type]
            if default_type in ["compound", "signature", "feature", "pathway", "dataset", "program", "experiment"]
            else 0,
        )
    with col2:
        entity_id_str = st.text_input("Entity UUID", value=default_id or "", placeholder="e.g. 2f1c... (UUID)")
    with col3:
        depth = st.slider("k-hop", min_value=1, max_value=3, value=2, step=1)

    relationships = st.multiselect(
        "Relationship filters",
        [
            "activity_against",
            "in_pathway",
            "enriched_in_pathway",
            "matches_signature",
        ],
        default=[],
        help="Optional: filter which edge types are traversed.",
    )

    col_opt1, col_opt2 = st.columns(2)
    with col_opt1:
        size_by_degree = st.checkbox("Size by Degree", value=False)
    with col_opt2:
        color_by_community = st.checkbox("Color by Community", value=False)

    if st.button("Traverse", type="primary"):
        try:
            entity_id = UUID(entity_id_str.strip())
        except Exception:
            st.error("Entity UUID invalid.")
            return

        payload: Dict[str, Any] = {"entity_type": entity_type, "entity_id": str(entity_id), "depth": depth}
        if relationships:
            payload["relationships"] = relationships

        try:
            with httpx.Client(timeout=15) as client:
                resp = client.post(TRAVERSE_ENDPOINT, json=payload)
            if resp.status_code >= 400:
                st.error(f"API error ({resp.status_code}): {resp.text}")
                return
            data = resp.json()
        except Exception as e:  # noqa: BLE001
            st.error(f"Request failed: {e}")
            return

        st.session_state["graph_explorer_subgraph"] = data

    data = st.session_state.get("graph_explorer_subgraph")
    if isinstance(data, dict) and data.get("cytoscape"):
        st.subheader("Neighborhood")
        if data.get("truncated"):
            st.warning("Result truncated to max nodes.")
        st.caption(f"Nodes: {len(data.get('nodes') or [])} | Edges: {len(data.get('edges') or [])}")

        degree = None
        comms = None
        if size_by_degree or color_by_community:
            metrics = []
            if size_by_degree:
                metrics.append("degree_centrality")
            if color_by_community:
                metrics.append("communities")
            try:
                with httpx.Client(timeout=15) as client:
                    resp = client.post(
                        ANALYTICS_ENDPOINT,
                        json={"nodes": data.get("nodes") or [], "edges": data.get("edges") or [], "metrics": metrics},
                    )
                if resp.status_code < 400:
                    a = resp.json() or {}
                    degree = a.get("degree_centrality")
                    comms = a.get("communities")
                    st.session_state["graph_explorer_analytics"] = a
                else:
                    st.warning("Analytics API returned error; using default styling.")
            except Exception:
                st.warning("Analytics unavailable; using default styling.")

        try:
            from streamlit_agraph import agraph, Config  # type: ignore[import-not-found]

            nodes, edges = _to_agraph(
                data["cytoscape"],
                degree_centrality=degree,
                communities=comms,
                size_by_degree=size_by_degree,
                color_by_community=color_by_community,
            )
            config = Config(
                width=900,
                height=600,
                directed=True,
                physics=True,
                hierarchical=False,
            )
            agraph(nodes=nodes, edges=edges, config=config)

            if color_by_community:
                analytics_data = st.session_state.get("graph_explorer_analytics", {})
                if analytics_data:
                    st.caption("**Community Detection:**")
                    for comm_id, size in (analytics_data.get("community_sizes") or {}).items():
                        st.caption(f"  - Community {comm_id}: {size} nodes")
                    mod = analytics_data.get("modularity")
                    if mod is not None:
                        st.metric("Modularity", f"{mod:.3f}", help="Higher = better community structure")
        except Exception as e:  # noqa: BLE001
            st.error(f"Visualization failed (is streamlit-agraph installed?): {e}")
            st.json(data["cytoscape"])

    st.divider()
    st.subheader("Shortest path")
    col_a, col_b = st.columns(2)
    with col_a:
        src_type = st.selectbox("Source type", ["compound", "signature", "feature", "pathway"], key="graph_src_type")
        src_id = st.text_input("Source UUID", key="graph_src_id")
    with col_b:
        tgt_type = st.selectbox("Target type", ["compound", "signature", "feature", "pathway"], key="graph_tgt_type")
        tgt_id = st.text_input("Target UUID", key="graph_tgt_id")

    rel2 = st.multiselect(
        "Relationship filters (path)",
        ["activity_against", "in_pathway", "enriched_in_pathway", "matches_signature"],
        default=[],
        key="graph_path_rels",
    )

    if st.button("Find path"):
        try:
            payload = {
                "source_type": src_type,
                "source_id": str(UUID(src_id.strip())),
                "target_type": tgt_type,
                "target_id": str(UUID(tgt_id.strip())),
            }
        except Exception:
            st.error("Source/target UUID invalid.")
            return
        if rel2:
            payload["relationships"] = rel2

        try:
            with httpx.Client(timeout=10) as client:
                resp = client.post(PATH_ENDPOINT, json=payload)
            if resp.status_code >= 400:
                st.error(f"API error ({resp.status_code}): {resp.text}")
                return
            out = resp.json()
        except Exception as e:  # noqa: BLE001
            st.error(f"Request failed: {e}")
            return

        if not out.get("found"):
            st.info("No path found (or timed out).")
        else:
            st.success(f"Path length: {max(0, len(out.get('nodes') or []) - 1)}")
            st.json(out.get("nodes"))


__all__ = ["render_graph_explorer_page"]


