from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import networkx as nx
import plotly.graph_objects as go
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Signature
from scripts.dashboard.db_session import db_session


def _list_signatures(db: Session) -> List[Dict[str, str]]:
    sigs = db.query(Signature).order_by(Signature.updated_at.desc()).limit(200).all()
    return [
        {"id": str(s.id), "name": str(s.name) if s.name is not None else str(s.id)}
        for s in sigs
        if s.id is not None
    ]


def _build_network(db: Session, signature_id: str) -> Tuple[List[str], List[Tuple[str, str]]]:
    sig = db.query(Signature).filter(Signature.id == signature_id).first()
    if not sig or not sig.features:
        return [], []
    nodes = [str(f.name) for f in sig.features if getattr(f, "name", None)]
    edges = list(itertools.combinations(nodes, 2))
    return nodes, edges


def _annotate_nodes_with_pathways(nodes: List[str]) -> Dict[str, str]:
    """Map feature names to their primary pathway (placeholder heuristic)."""
    pathways: Dict[str, str] = {}
    for node in nodes:
        if "TP53" in node or "MDM2" in node:
            pathways[node] = "Cell Cycle"
        elif "AKT" in node or "PI3K" in node:
            pathways[node] = "PI3K-Akt Signaling"
        else:
            pathways[node] = "Other"
    return pathways


def render() -> None:
    st.header("Signature Network")
    st.caption("Feature co-occurrence network from Postgres signatures.")

    with db_session() as db:
        signatures = _list_signatures(db)

    if not signatures:
        st.info("No signatures available.")
        return

    sig_options = {s["name"]: s["id"] for s in signatures if s.get("id")}
    sig_label = st.selectbox("Signature", list(sig_options.keys()), index=0, key="network_signature")
    sig_id = sig_options[sig_label]

    if st.button("Refresh data", key="network_refresh"):
        st.rerun()

    if st.button("Build Network", key="network_generate"):
        with st.spinner("Building network from Postgres..."):
            with db_session() as db:
                nodes, edges = _build_network(db, sig_id)

        if not nodes:
            st.warning("No features found for this signature.")
            return

        pathway_map = _annotate_nodes_with_pathways(nodes)
        unique_paths = sorted(set(pathway_map.values()))
        pathway_filter = st.selectbox(
            "Filter by Pathway",
            ["All"] + unique_paths,
            key="network_pathway_filter",
        )

        if pathway_filter != "All":
            filtered_nodes = [n for n in nodes if pathway_map.get(n, "Other") == pathway_filter]
            filtered_edges = [(s, t) for (s, t) in edges if s in filtered_nodes and t in filtered_nodes]
        else:
            filtered_nodes = nodes
            filtered_edges = edges

        st.metric("Nodes", len(filtered_nodes))
        st.metric("Edges", len(filtered_edges))

        # Build force-directed layout
        G = nx.Graph()
        G.add_nodes_from(filtered_nodes)
        G.add_edges_from(filtered_edges)
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

        color_palette = {
            "Cell Cycle": "#E15759",
            "PI3K-Akt Signaling": "#4E79A7",
            "Other": "#59A14F",
        }
        colors = [
            color_palette.get(pathway_map.get(n, "Other"), "#59A14F") for n in filtered_nodes
        ]

        edge_x, edge_y = [], []
        for s, t in filtered_edges:
            x0, y0 = pos[s]
            x1, y1 = pos[t]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]

        node_x = [pos[n][0] for n in filtered_nodes]
        node_y = [pos[n][1] for n in filtered_nodes]

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=edge_x,
                y=edge_y,
                mode="lines",
                line=dict(color="#888", width=1),
                hoverinfo="none",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                text=filtered_nodes,
                textposition="top center",
                marker=dict(size=16, color=colors),
            )
        )
        fig.update_layout(showlegend=False)

        st.plotly_chart(fig, width='stretch')

        st.caption(
            "Legend: 🔴 Cell Cycle  |  🔵 PI3K-Akt Signaling  |  ⚫ Other"
        )

