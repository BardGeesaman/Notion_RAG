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


def _to_agraph(cytoscape: Dict[str, Any]):
    from streamlit_agraph import Node, Edge  # type: ignore[import-not-found]

    COLORS = {
        "compound": "#2E86AB",
        "signature": "#A23B72",
        "feature": "#F18F01",
        "pathway": "#C73E1D",
        "dataset": "#3B1F2B",
    }

    nodes = []
    for n in cytoscape.get("nodes", []):
        d = (n or {}).get("data") or {}
        color = COLORS.get(d.get("entity_type"), "#888888")
        nodes.append(
            Node(
                id=d.get("id"),
                label=d.get("label"),
                size=20,
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

        try:
            from streamlit_agraph import agraph, Config  # type: ignore[import-not-found]

            nodes, edges = _to_agraph(data["cytoscape"])
            config = Config(
                width=900,
                height=600,
                directed=True,
                physics=True,
                hierarchical=False,
            )
            agraph(nodes=nodes, edges=edges, config=config)
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


