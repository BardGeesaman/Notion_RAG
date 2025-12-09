"""Relationship visualization page for the Streamlit dashboard."""

from __future__ import annotations

import uuid
from typing import Dict, Set, Tuple

import streamlit as st

try:
    import networkx as nx
    import plotly.graph_objects as go

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

from amprenta_rag.database.models import Dataset, Experiment, Feature, Program, Signature
from scripts.dashboard.db_session import db_session


def render_relationships_page() -> None:
    """
    Render the Relationship Visualization page.

    Displays interactive network graphs showing relationships between:
    - Programs, Datasets, Experiments
    - Datasets, Features
    - Programs, Signatures
    - Datasets, Signatures
    """
    st.header("🔗 Relationship Visualization")

    if not HAS_NETWORKX:
        st.error("⚠️ Network visualization requires networkx and plotly. " "Install with: `pip install networkx plotly`")
        return

    with db_session() as db:
        # Entity selection
        st.subheader("Select Entities to Visualize")

        col1, col2 = st.columns(2)

        with col1:
            entity_type = st.selectbox(
                "Start with:",
                ["Program", "Dataset", "Experiment", "Feature", "Signature"],
                help="Select the type of entity to start the visualization from",
            )

        with col2:
            max_depth = st.slider(
                "Relationship Depth:",
                min_value=1,
                max_value=3,
                value=2,
                help="How many relationship levels to show",
            )

        # Get all entities for selection
        if entity_type == "Program":
            entities = db.query(Program).all()
            entity_list = [(str(p.id), p.name) for p in entities]
        elif entity_type == "Dataset":
            entities = db.query(Dataset).all()
            entity_list = [(str(d.id), f"{d.name} ({d.omics_type})") for d in entities]
        elif entity_type == "Experiment":
            entities = db.query(Experiment).all()
            entity_list = [(str(e.id), e.name) for e in entities]
        elif entity_type == "Feature":
            entities = db.query(Feature).limit(100).all()  # Limit features
            entity_list = [(str(f.id), f"{f.name} ({f.feature_type})") for f in entities]
        else:  # Signature
            entities = db.query(Signature).all()
            entity_list = [(str(s.id), s.name) for s in entities]

        if not entity_list:
            st.info(f"No {entity_type.lower()}s found in the database.")
            return

        selected_entity_id = st.selectbox(
            f"Select {entity_type}:",
            options=[e[0] for e in entity_list],
            format_func=lambda x: next((e[1] for e in entity_list if e[0] == x), x),
        )

        if st.button("🔍 Visualize Relationships", type="primary"):
            with st.spinner("Building relationship graph..."):
                graph, node_info = build_relationship_graph(
                    db=db,
                    start_type=entity_type,
                    start_id=uuid.UUID(selected_entity_id),
                    max_depth=max_depth,
                )

                if graph.number_of_nodes() == 0:
                    st.warning("No relationships found for the selected entity.")
                    return

                # Display graph
                st.subheader("📊 Relationship Graph")
                fig = create_network_visualization(graph, node_info)
                st.plotly_chart(fig, width='stretch')

                # Display statistics
                st.subheader("📈 Graph Statistics")
                col1, col2, col3, col4 = st.columns(4)

                with col1:
                    st.metric("Nodes", graph.number_of_nodes())
                with col2:
                    st.metric("Edges", graph.number_of_edges())
                with col3:
                    st.metric("Components", nx.number_connected_components(graph))
                with col4:
                    if graph.number_of_nodes() > 0:
                        density = nx.density(graph)
                        st.metric("Density", f"{density:.3f}")

                # Display node details
                st.subheader("🔍 Node Details")
                node_types = {}
                for node_id, info in node_info.items():
                    node_type = info["type"]
                    if node_type not in node_types:
                        node_types[node_type] = []
                    node_types[node_type].append(info["name"])

                for node_type, names in sorted(node_types.items()):
                    with st.expander(f"{node_type}s ({len(names)})"):
                        for name in sorted(names)[:50]:  # Limit display
                            st.write(f"• {name}")
                        if len(names) > 50:
                            st.write(f"... and {len(names) - 50} more")


def build_relationship_graph(
    db,
    start_type: str,
    start_id: uuid.UUID,
    max_depth: int,
) -> Tuple[nx.Graph, Dict[str, Dict]]:
    """
    Build a NetworkX graph from database relationships.

    Args:
        db: Database session
        start_type: Type of starting entity (Program, Dataset, etc.)
        start_id: UUID of starting entity
        max_depth: Maximum relationship depth to traverse

    Returns:
        Tuple of (NetworkX graph, node_info dictionary)
    """
    graph = nx.Graph()
    node_info: Dict[str, Dict] = {}
    visited: Set[Tuple[str, uuid.UUID]] = set()

    def add_node(entity_type: str, entity_id: uuid.UUID, name: str, **kwargs):
        """Add a node to the graph."""
        node_id = f"{entity_type}:{entity_id}"
        if node_id not in graph:
            graph.add_node(node_id)
            node_info[node_id] = {
                "type": entity_type,
                "id": str(entity_id),
                "name": name,
                **kwargs,
            }
        return node_id

    def add_edge(node1: str, node2: str, edge_type: str):
        """Add an edge to the graph."""
        if not graph.has_edge(node1, node2):
            graph.add_edge(node1, node2, type=edge_type)

    def traverse_relationships(
        entity_type: str,
        entity_id: uuid.UUID,
        depth: int,
    ):
        """Recursively traverse relationships."""
        if depth > max_depth:
            return

        visit_key = (entity_type, entity_id)
        if visit_key in visited:
            return
        visited.add(visit_key)

        if entity_type == "Program":
            program = db.query(Program).filter(Program.id == entity_id).first()
            if not program:
                return

            node_id = add_node("Program", program.id, program.name)

            # Programs → Datasets
            for dataset in program.datasets:
                ds_node = add_node("Dataset", dataset.id, dataset.name, omics_type=dataset.omics_type)
                add_edge(node_id, ds_node, "has_dataset")
                if depth < max_depth:
                    traverse_relationships("Dataset", dataset.id, depth + 1)

            # Programs → Experiments
            for experiment in program.experiments:
                exp_node = add_node("Experiment", experiment.id, experiment.name)
                add_edge(node_id, exp_node, "has_experiment")
                if depth < max_depth:
                    traverse_relationships("Experiment", experiment.id, depth + 1)

            # Programs → Signatures
            for signature in program.signatures:
                sig_node = add_node("Signature", signature.id, signature.name)
                add_edge(node_id, sig_node, "has_signature")

        elif entity_type == "Dataset":
            dataset = db.query(Dataset).filter(Dataset.id == entity_id).first()
            if not dataset:
                return

            node_id = add_node("Dataset", dataset.id, dataset.name, omics_type=dataset.omics_type)

            # Datasets → Programs
            for program in dataset.programs:
                prog_node = add_node("Program", program.id, program.name)
                add_edge(node_id, prog_node, "belongs_to_program")
                if depth < max_depth:
                    traverse_relationships("Program", program.id, depth + 1)

            # Datasets → Experiments
            for experiment in dataset.experiments:
                exp_node = add_node("Experiment", experiment.id, experiment.name)
                add_edge(node_id, exp_node, "belongs_to_experiment")
                if depth < max_depth:
                    traverse_relationships("Experiment", experiment.id, depth + 1)

            # Datasets → Features
            for feature in dataset.features:
                feat_node = add_node("Feature", feature.id, feature.name, feature_type=feature.feature_type)
                add_edge(node_id, feat_node, "has_feature")

            # Datasets → Signatures
            for signature in dataset.signatures:
                sig_node = add_node("Signature", signature.id, signature.name)
                add_edge(node_id, sig_node, "matches_signature")

        elif entity_type == "Experiment":
            experiment = db.query(Experiment).filter(Experiment.id == entity_id).first()
            if not experiment:
                return

            node_id = add_node("Experiment", experiment.id, experiment.name)

            # Experiments → Programs
            for program in experiment.programs:
                prog_node = add_node("Program", program.id, program.name)
                add_edge(node_id, prog_node, "belongs_to_program")
                if depth < max_depth:
                    traverse_relationships("Program", program.id, depth + 1)

            # Experiments → Datasets
            for dataset in experiment.datasets:
                ds_node = add_node("Dataset", dataset.id, dataset.name, omics_type=dataset.omics_type)
                add_edge(node_id, ds_node, "has_dataset")
                if depth < max_depth:
                    traverse_relationships("Dataset", dataset.id, depth + 1)

        elif entity_type == "Feature":
            feature = db.query(Feature).filter(Feature.id == entity_id).first()
            if not feature:
                return

            node_id = add_node("Feature", feature.id, feature.name, feature_type=feature.feature_type)

            # Features → Datasets
            for dataset in feature.datasets:
                ds_node = add_node("Dataset", dataset.id, dataset.name, omics_type=dataset.omics_type)
                add_edge(node_id, ds_node, "in_dataset")
                if depth < max_depth:
                    traverse_relationships("Dataset", dataset.id, depth + 1)

            # Features → Signatures
            for signature in feature.signatures:
                sig_node = add_node("Signature", signature.id, signature.name)
                add_edge(node_id, sig_node, "in_signature")

        elif entity_type == "Signature":
            signature = db.query(Signature).filter(Signature.id == entity_id).first()
            if not signature:
                return

            node_id = add_node("Signature", signature.id, signature.name)

            # Signatures → Programs
            for program in signature.programs:
                prog_node = add_node("Program", program.id, program.name)
                add_edge(node_id, prog_node, "in_program")
                if depth < max_depth:
                    traverse_relationships("Program", program.id, depth + 1)

            # Signatures → Datasets
            for dataset in signature.datasets:
                ds_node = add_node("Dataset", dataset.id, dataset.name, omics_type=dataset.omics_type)
                add_edge(node_id, ds_node, "matched_by_dataset")
                if depth < max_depth:
                    traverse_relationships("Dataset", dataset.id, depth + 1)

            # Signatures → Features
            for feature in signature.features:
                feat_node = add_node("Feature", feature.id, feature.name, feature_type=feature.feature_type)
                add_edge(node_id, feat_node, "has_feature")

    # Start traversal
    traverse_relationships(start_type, start_id, depth=0)

    return graph, node_info


def create_network_visualization(
    graph: nx.Graph,
    node_info: Dict[str, Dict],
) -> go.Figure:
    """
    Create an interactive Plotly network visualization.

    Args:
        graph: NetworkX graph
        node_info: Dictionary mapping node IDs to node information

    Returns:
        Plotly figure object
    """
    # Use spring layout for positioning
    pos = nx.spring_layout(graph, k=1, iterations=50)

    # Color mapping for node types
    type_colors = {
        "Program": "#FF6B6B",
        "Dataset": "#4ECDC4",
        "Experiment": "#95E1D3",
        "Feature": "#F38181",
        "Signature": "#AA96DA",
    }

    # Create edge traces
    edge_traces = []
    for edge in graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_traces.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode="lines",
                line=dict(width=1, color="#888"),
                hoverinfo="none",
                showlegend=False,
            )
        )

    # Create node traces by type
    node_traces = []
    for node_type in ["Program", "Dataset", "Experiment", "Feature", "Signature"]:
        nodes_of_type = [node for node in graph.nodes() if node_info[node]["type"] == node_type]

        if not nodes_of_type:
            continue

        x_coords = [pos[node][0] for node in nodes_of_type]
        y_coords = [pos[node][1] for node in nodes_of_type]
        node_names = [node_info[node]["name"] for node in nodes_of_type]
        node_ids = [node_info[node]["id"] for node in nodes_of_type]

        node_traces.append(
            go.Scatter(
                x=x_coords,
                y=y_coords,
                mode="markers+text",
                name=node_type,
                marker=dict(
                    size=15,
                    color=type_colors.get(node_type, "#999"),
                    line=dict(width=2, color="white"),
                ),
                text=[name[:30] + "..." if len(name) > 30 else name for name in node_names],
                textposition="middle center",
                textfont=dict(size=8),
                hovertemplate="<b>%{text}</b><br>Type: %{fullData.name}<br>ID: %{customdata}<extra></extra>",
                customdata=node_ids,
            )
        )

    # Combine all traces
    all_traces = edge_traces + node_traces

    # Create figure
    fig = go.Figure(
        data=all_traces,
        layout=go.Layout(
            title="",
            showlegend=True,
            hovermode="closest",
            margin=dict(b=20, l=5, r=5, t=40),
            annotations=[
                dict(
                    text="Interactive network graph - hover over nodes for details",
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.005,
                    y=-0.002,
                    xanchor="left",
                    yanchor="bottom",
                    font=dict(size=12, color="#888"),
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
        ),
    )

    return fig
