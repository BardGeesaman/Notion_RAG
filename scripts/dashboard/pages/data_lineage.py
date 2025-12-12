"""Data Lineage Visualization page."""
from __future__ import annotations

from uuid import UUID

import plotly.graph_objects as go
import streamlit as st

from amprenta_rag.database.models import Experiment, Dataset, Signature
from amprenta_rag.utils.lineage import get_entity_lineage, get_full_lineage_graph
from scripts.dashboard.db_session import db_session


def render_data_lineage_page() -> None:
    """Render the Data Lineage Visualization page."""
    st.header("🔗 Data Lineage")
    st.markdown("Visualize relationships between experiments, datasets, signatures, and protocols")
    
    # Mode selector
    view_mode = st.radio("View", ["Single Entity", "Full Graph"], horizontal=True)
    
    with db_session() as db:
        if view_mode == "Single Entity":
            render_single_entity_view(db)
        else:
            render_full_graph_view(db)


def render_single_entity_view(db) -> None:
    """Render single entity lineage view."""
    st.subheader("Single Entity Lineage")
    
    # Entity type selector
    entity_type = st.selectbox("Entity Type", ["Experiment", "Dataset", "Signature"])
    
    # Entity selector based on type
    entity_id = None
    entity_name = None
    
    if entity_type == "Experiment":
        experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()
        if experiments:
            experiment_options = {exp.name: exp.id for exp in experiments}
            selected_name = st.selectbox("Select Experiment", list(experiment_options.keys()))
            entity_id = experiment_options[selected_name]
            entity_name = selected_name
        else:
            st.info("No experiments available.")
            return
    elif entity_type == "Dataset":
        datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()
        if datasets:
            dataset_options = {ds.name: ds.id for ds in datasets}
            selected_name = st.selectbox("Select Dataset", list(dataset_options.keys()))
            entity_id = dataset_options[selected_name]
            entity_name = selected_name
        else:
            st.info("No datasets available.")
            return
    elif entity_type == "Signature":
        signatures = db.query(Signature).order_by(Signature.name).limit(200).all()
        if signatures:
            signature_options = {sig.name: sig.id for sig in signatures}
            selected_name = st.selectbox("Select Signature", list(signature_options.keys()))
            entity_id = signature_options[selected_name]
            entity_name = selected_name
        else:
            st.info("No signatures available.")
            return
    
    if st.button("Show Lineage", type="primary"):
        if entity_id:
            with st.spinner("Building lineage graph..."):
                try:
                    lineage_data = get_entity_lineage(entity_type.lower(), entity_id, db)
                    render_lineage_graph(lineage_data, title=f"Lineage: {entity_name}")
                except Exception as e:
                    st.error(f"Error building lineage: {e}")
                    st.exception(e)


def render_full_graph_view(db) -> None:
    """Render full graph overview."""
    st.subheader("Full Lineage Graph")
    st.markdown("Overview of all connections in the system")
    
    limit = st.slider("Max Experiments", min_value=10, max_value=500, value=100, step=10)
    
    if st.button("Generate Graph", type="primary"):
        with st.spinner(f"Building full lineage graph (up to {limit} experiments)..."):
            try:
                lineage_data = get_full_lineage_graph(db, limit=limit)
                render_lineage_graph(lineage_data, title="Full Lineage Graph")
            except Exception as e:
                st.error(f"Error building full graph: {e}")
                st.exception(e)


def render_lineage_graph(lineage_data: dict, title: str = "Lineage Graph") -> None:
    """
    Render a lineage graph using Plotly.
    
    Args:
        lineage_data: Dictionary with "nodes" and "edges" lists
        title: Graph title
    """
    nodes = lineage_data.get("nodes", [])
    edges = lineage_data.get("edges", [])
    
    if not nodes:
        st.info("No nodes found in lineage graph.")
        return
    
    if not edges:
        st.warning("No connections found. Graph may show isolated nodes.")
    
    # Color mapping for node types
    color_map = {
        "experiment": "#3498db",  # Blue
        "dataset": "#2ecc71",      # Green
        "signature": "#e67e22",    # Orange
        "protocol": "#9b59b6",     # Purple
        "project": "#f39c12",      # Yellow
    }
    
    # Build node positions using a simple force-directed layout simulation
    # For simplicity, we'll use a circular layout
    import math
    
    node_positions = {}
    num_nodes = len(nodes)
    
    # Calculate positions in a circle
    for i, node in enumerate(nodes):
        angle = 2 * math.pi * i / num_nodes if num_nodes > 0 else 0
        radius = 3.0
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        node_positions[node["id"]] = (x, y)
    
    # Create edge traces
    edge_x = []
    edge_y = []
    edge_info = []
    
    for edge in edges:
        source_pos = node_positions.get(edge["source"])
        target_pos = node_positions.get(edge["target"])
        
        if source_pos and target_pos:
            edge_x.extend([source_pos[0], target_pos[0], None])
            edge_y.extend([source_pos[1], target_pos[1], None])
            edge_info.append(edge["relationship"])
    
    # Create node traces
    node_x = []
    node_y = []
    node_text = []
    node_colors = []
    node_hover = []
    
    for node in nodes:
        pos = node_positions.get(node["id"])
        if pos:
            node_x.append(pos[0])
            node_y.append(pos[1])
            node_text.append(node["name"])
            node_colors.append(color_map.get(node["type"], "#95a5a6"))  # Default gray
            
            # Build hover text
            hover_parts = [f"<b>{node['name']}</b>", f"Type: {node['type']}"]
            if node.get("metadata"):
                for key, value in node["metadata"].items():
                    if value:
                        hover_parts.append(f"{key}: {value}")
            node_hover.append("<br>".join(hover_parts))
    
    # Create figure
    fig = go.Figure()
    
    # Add edges
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(width=1, color="#888"),
            hoverinfo="none",
            showlegend=False,
        )
    )
    
    # Add nodes
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            marker=dict(
                size=20,
                color=node_colors,
                line=dict(width=2, color="white"),
            ),
            text=node_text,
            textposition="middle center",
            textfont=dict(size=10, color="white"),
            hovertext=node_hover,
            hoverinfo="text",
            showlegend=False,
        )
    )
    
    # Update layout
    fig.update_layout(
        title=title,
        showlegend=False,
        hovermode="closest",
        margin=dict(b=20, l=5, r=5, t=40),
        annotations=[
            dict(
                text="Node colors: 🔵 Experiment, 🟢 Dataset, 🟠 Signature, 🟣 Protocol, 🟡 Project",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.005,
                y=-0.002,
                xanchor="left",
                yanchor="bottom",
                font=dict(size=10, color="#666"),
            )
        ],
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600,
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Show summary stats
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Nodes", len(nodes))
    with col2:
        st.metric("Total Edges", len(edges))
    with col3:
        node_types = {}
        for node in nodes:
            node_types[node["type"]] = node_types.get(node["type"], 0) + 1
        st.metric("Node Types", len(node_types))
    
    # Show node type breakdown
    if node_types:
        st.markdown("### Node Type Breakdown")
        type_data = {k: v for k, v in sorted(node_types.items(), key=lambda x: x[1], reverse=True)}
        st.bar_chart(type_data)

