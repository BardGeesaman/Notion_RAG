from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import plotly.graph_objects as go
import streamlit as st
from sqlalchemy.orm import Session

from amprenta_rag.database.models import Signature
from scripts.dashboard.db_session import db_session


def _list_signatures(db: Session) -> List[Dict[str, str]]:
    sigs = db.query(Signature).order_by(Signature.updated_at.desc()).limit(200).all()
    return [{"id": str(s.id), "name": s.name} for s in sigs]


def _build_network(db: Session, signature_id: str) -> Tuple[List[str], List[Tuple[str, str]]]:
    sig = db.query(Signature).filter(Signature.id == signature_id).first()
    if not sig or not sig.features:
        return [], []
    nodes = [f.name for f in sig.features]
    edges = list(itertools.combinations(nodes, 2))
    return nodes, edges


def render() -> None:
    st.header("Signature Network")
    st.caption("Feature co-occurrence network from Postgres signatures.")

    with db_session() as db:
        signatures = _list_signatures(db)

    if not signatures:
        st.info("No signatures available.")
        return

    sig_options = {s["name"]: s["id"] for s in signatures}
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

        # simple grid layout
        pos = {n: (i % 6, i // 6) for i, n in enumerate(nodes)}
        edge_x, edge_y = [], []
        for s, t in edges:
            x0, y0 = pos[s]
            x1, y1 = pos[t]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]

        node_x = [pos[n][0] for n in nodes]
        node_y = [pos[n][1] for n in nodes]

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
                text=nodes,
                textposition="top center",
                marker=dict(size=16, color="#4E79A7"),
            )
        )
        fig.update_layout(showlegend=False)

        st.plotly_chart(fig, width='stretch')

