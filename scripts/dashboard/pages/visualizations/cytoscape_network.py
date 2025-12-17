from __future__ import annotations

import streamlit as st

from scripts.dashboard.components.cytoscape import render_cytoscape


def render() -> None:
    st.header("Cytoscape Network")
    st.caption("Interactive Cytoscape.js network visualization.")

    st.markdown("### Sample Network")
    nodes = [
        {"id": "A", "label": "Node A"},
        {"id": "B", "label": "Node B"},
        {"id": "C", "label": "Node C"},
        {"id": "D", "label": "Node D"},
    ]
    edges = [
        {"source": "A", "target": "B"},
        {"source": "A", "target": "C"},
        {"source": "B", "target": "D"},
        {"source": "C", "target": "D"},
    ]

    layout = st.selectbox("Layout", ["cose", "grid", "circle", "breadthfirst"], index=0)
    height = st.slider("Height (px)", min_value=300, max_value=900, value=600, step=50)

    render_cytoscape(nodes, edges, height=height, layout=layout)

    st.info(
        "Click nodes or edges to see selection below the graph. "
        "Pass your own nodes/edges to render_cytoscape for custom networks."
    )

