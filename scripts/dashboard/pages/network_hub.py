"""Unified Network Hub for Cytoscape.js visualizations."""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

import httpx
import streamlit as st

from scripts.dashboard.components.cytoscape import render_cytoscape


def _api_base() -> str:
    return os.environ.get("API_URL", "http://localhost:8000").rstrip("/")


def _api_get(path: str, params: Optional[Dict] = None) -> Any:
    with httpx.Client(timeout=60) as client:
        r = client.get(f"{_api_base()}{path}", params=params)
        r.raise_for_status()
        return r.json()


def _api_post(path: str, payload: Dict) -> Any:
    with httpx.Client(timeout=60) as client:
        r = client.post(f"{_api_base()}{path}", json=payload)
        r.raise_for_status()
        return r.json()


@st.cache_data(ttl=300)
def _list_datasets(limit: int = 100) -> List[Dict]:
    try:
        return _api_get("/api/v1/datasets", {"limit": limit})
    except Exception:
        return []


def render_network_hub_page() -> None:
    """Main Network Hub page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🕸️ Network Hub")
    st.caption("Unified network visualization: PPI, Evidence Graph, Compound-Target, Pathways")
    
    # Shared controls in sidebar
    with st.sidebar:
        st.subheader("Network Controls")
        layout = st.selectbox("Layout", ["cose", "concentric", "circle", "grid"], index=0)
        detect_clusters = st.checkbox("Detect Clusters", value=False)
        show_edge_labels = st.checkbox("Show Edge Labels", value=False)
    
    # Tabs
    tab_ppi, tab_evidence, tab_compound, tab_pathways = st.tabs([
        "🧬 PPI Networks",
        "🕸️ Evidence Graph",
        "💊 Compound-Target",
        "🗺️ Pathways"
    ])
    
    with tab_ppi:
        _render_ppi_tab(layout, detect_clusters, show_edge_labels)
    
    with tab_evidence:
        _render_evidence_tab(layout, detect_clusters, show_edge_labels)
    
    with tab_compound:
        _render_compound_tab(layout, detect_clusters, show_edge_labels)
    
    with tab_pathways:
        _render_pathways_tab()


def _render_ppi_tab(layout: str, detect_clusters: bool, show_edge_labels: bool) -> None:
    """PPI Networks tab using STRING API."""
    st.subheader("Protein-Protein Interaction Networks")
    st.caption("Build networks from STRING database")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        genes_input = st.text_area(
            "Gene Symbols",
            placeholder="Enter gene symbols (one per line or comma-separated)\nExample:\nTP53\nBRCA1\nMDM2",
            height=150,
            key="ppi_genes"
        )
    
    with col2:
        species = st.selectbox(
            "Species",
            options=[9606, 10090],
            format_func=lambda x: "Human" if x == 9606 else "Mouse",
            key="ppi_species"
        )
        min_score = st.slider("Confidence", 0, 1000, 400, 50, key="ppi_score")
        st.caption(f"Score ≥ {min_score} ({['Low', 'Medium', 'High'][min(min_score // 350, 2)]} confidence)")
    
    # Parse genes
    genes = []
    if genes_input:
        genes = [g.strip().upper() for g in genes_input.replace(",", "\n").split("\n") if g.strip()]
    
    if st.button("🔨 Build PPI Network", disabled=len(genes) < 2, key="ppi_build"):
        with st.spinner(f"Fetching PPI network for {len(genes)} genes..."):
            try:
                result = _api_post("/api/graph/ppi/network", {
                    "proteins": genes,
                    "species": species,
                    "min_score": min_score,
                })
                st.session_state["ppi_network"] = result
                st.success(f"Found {result.get('interaction_count', 0)} interactions")
            except Exception as e:
                st.error(f"Failed to fetch PPI network: {e}")
    
    # Expression overlay section
    st.divider()
    st.subheader("Expression Overlay")
    
    datasets = _list_datasets()
    dataset_options = {d.get("name", str(d.get("id"))[:8]): d.get("id") for d in datasets if d.get("id")}
    
    col1, col2 = st.columns([2, 1])
    with col1:
        selected_dataset = st.selectbox(
            "Dataset",
            options=list(dataset_options.keys()) if dataset_options else ["No datasets available"],
            key="ppi_dataset"
        )
    with col2:
        colormap = st.selectbox("Colormap", ["diverging", "sequential"], key="ppi_colormap")
    
    apply_overlay = st.button("🎨 Apply Expression Overlay", key="ppi_overlay")
    
    # Render network
    network = st.session_state.get("ppi_network")
    if network:
        nodes = network.get("nodes", [])
        edges = network.get("edges", [])
        
        node_colors = None
        communities = None
        color_by = "default"
        
        # Apply expression overlay
        if apply_overlay and selected_dataset in dataset_options:
            dataset_id = dataset_options[selected_dataset]
            gene_symbols = [n.get("data", {}).get("id") or n.get("id") for n in nodes]
            try:
                overlay = _api_post(f"/api/v1/datasets/{dataset_id}/expression-overlay", {
                    "gene_symbols": gene_symbols,
                    "colormap": colormap,
                })
                node_colors = overlay.get("node_colors", {})
                color_by = "expression"
                st.info(f"Expression data: {overlay.get('genes_found', 0)} genes found")
            except Exception as e:
                st.warning(f"Expression overlay failed: {e}")
        
        # Detect clusters
        if detect_clusters and nodes:
            try:
                analytics = _api_post("/api/graph/analytics", {
                    "nodes": nodes,
                    "edges": edges,
                    "metrics": ["communities"],
                })
                communities = analytics.get("communities", {})
                if communities and color_by == "default":
                    color_by = "community"
            except Exception:
                pass
        
        st.metric("Nodes", len(nodes))
        st.metric("Edges", len(edges))
        
        render_cytoscape(
            nodes=nodes,
            edges=edges,
            height=600,
            layout=layout,
            node_colors=node_colors,
            communities=communities,
            show_edge_labels=show_edge_labels,
            color_by=color_by,
        )


def _render_evidence_tab(layout: str, detect_clusters: bool, show_edge_labels: bool) -> None:
    """Evidence Graph tab - full Graph Explorer functionality."""
    st.subheader("Evidence Graph Explorer")
    st.caption("Explore entity neighborhoods and shortest paths")
    
    # Deep-link support
    qp = getattr(st, "query_params", {})
    default_type = qp.get("entity_type", [None])[0] if isinstance(qp.get("entity_type"), list) else qp.get("entity_type")
    default_id = qp.get("entity_id", [None])[0] if isinstance(qp.get("entity_id"), list) else qp.get("entity_id")
    
    # Mode selection
    mode = st.radio("Mode", ["Traverse", "Shortest Path"], horizontal=True, key="evidence_mode")
    
    if mode == "Traverse":
        _render_traverse_mode(layout, detect_clusters, show_edge_labels, default_type, default_id)
    else:
        _render_path_mode()


def _render_traverse_mode(layout: str, detect_clusters: bool, show_edge_labels: bool, 
                          default_type: str = None, default_id: str = None) -> None:
    """K-hop neighborhood traversal."""
    ENTITY_TYPES = ["compound", "signature", "feature", "pathway", "dataset", "program", "experiment"]
    
    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        default_idx = ENTITY_TYPES.index(default_type) if default_type in ENTITY_TYPES else 0
        entity_type = st.selectbox("Entity type", ENTITY_TYPES, index=default_idx, key="ev_entity_type")
    with col2:
        entity_id_str = st.text_input("Entity UUID", value=default_id or "", key="ev_entity_id")
    with col3:
        depth = st.slider("k-hop depth", 1, 3, 2, key="ev_depth")
    
    # Relationship filters
    relationships = st.multiselect(
        "Relationship filters",
        ["activity_against", "in_pathway", "enriched_in_pathway", "matches_signature"],
        default=[],
        key="ev_relationships"
    )
    
    # Visualization options
    col_opt1, col_opt2 = st.columns(2)
    with col_opt1:
        size_by_degree = st.checkbox("Size by Degree", value=False, key="ev_size_degree")
    with col_opt2:
        # Use sidebar's detect_clusters or local toggle
        color_by_community = detect_clusters or st.checkbox("Color by Community", value=False, key="ev_color_comm")
    
    if st.button("🔍 Traverse", type="primary", key="ev_traverse_btn"):
        try:
            from uuid import UUID
            entity_id = UUID(entity_id_str.strip())
        except Exception:
            st.error("Invalid UUID format")
            return
        
        payload = {
            "entity_type": entity_type,
            "entity_id": str(entity_id),
            "depth": depth,
        }
        if relationships:
            payload["relationships"] = relationships
        
        with st.spinner("Traversing graph..."):
            try:
                result = _api_post("/api/graph/traverse", payload)
                st.session_state["evidence_subgraph"] = result
                st.success(f"Found {len(result.get('nodes', []))} nodes, {len(result.get('edges', []))} edges")
            except Exception as e:
                st.error(f"Traversal failed: {e}")
    
    # Render network
    data = st.session_state.get("evidence_subgraph")
    if data and data.get("cytoscape"):
        if data.get("truncated"):
            st.warning("Result truncated to max nodes")
        
        cytoscape = data["cytoscape"]
        nodes = cytoscape.get("nodes", [])
        edges = cytoscape.get("edges", [])
        
        # Analytics for sizing/coloring
        degree = None
        communities = None
        
        if size_by_degree or color_by_community:
            metrics = []
            if size_by_degree:
                metrics.append("degree_centrality")
            if color_by_community:
                metrics.append("communities")
            
            try:
                analytics = _api_post("/api/graph/analytics", {
                    "nodes": data.get("nodes", []),
                    "edges": data.get("edges", []),
                    "metrics": metrics,
                })
                degree = analytics.get("degree_centrality")
                communities = analytics.get("communities")
                st.session_state["evidence_analytics"] = analytics
            except Exception:
                st.warning("Analytics unavailable")
        
        # Apply degree sizing to nodes
        if size_by_degree and degree:
            max_deg = max(degree.values()) if degree else 1.0
            for n in nodes:
                nid = n.get("data", {}).get("id")
                if nid and nid in degree:
                    # Scale size: 20-60px based on degree
                    n["data"]["size"] = 20 + int(40 * (degree[nid] / max_deg))
        
        # Determine color mode
        color_by = "default"
        if color_by_community and communities:
            color_by = "community"
        
        # Show metrics
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Nodes", len(nodes))
        with col2:
            st.metric("Edges", len(edges))
        
        # Show community stats if available
        if color_by_community:
            analytics = st.session_state.get("evidence_analytics", {})
            if analytics.get("community_sizes"):
                with st.expander("Community Statistics"):
                    for cid, size in analytics["community_sizes"].items():
                        st.caption(f"Community {cid}: {size} nodes")
                    if analytics.get("modularity"):
                        st.metric("Modularity", f"{analytics['modularity']:.3f}")
        
        render_cytoscape(
            nodes=nodes,
            edges=edges,
            height=550,
            layout=layout,
            communities=communities,
            show_edge_labels=show_edge_labels,
            color_by=color_by,
        )
        
        # Export options
        with st.expander("Export"):
            st.download_button(
                "📥 Download JSON",
                data=json.dumps(cytoscape, indent=2),
                file_name="evidence_graph.json",
                mime="application/json",
            )


def _render_path_mode() -> None:
    """Shortest path finding."""
    st.subheader("Shortest Path")
    
    ENTITY_TYPES = ["compound", "signature", "feature", "pathway"]
    
    col_a, col_b = st.columns(2)
    with col_a:
        src_type = st.selectbox("Source type", ENTITY_TYPES, key="ev_src_type")
        src_id = st.text_input("Source UUID", key="ev_src_id")
    with col_b:
        tgt_type = st.selectbox("Target type", ENTITY_TYPES, key="ev_tgt_type")
        tgt_id = st.text_input("Target UUID", key="ev_tgt_id")
    
    rel_filter = st.multiselect(
        "Relationship filters",
        ["activity_against", "in_pathway", "enriched_in_pathway", "matches_signature"],
        default=[],
        key="ev_path_rels"
    )
    
    if st.button("🔎 Find Path", key="ev_path_btn"):
        try:
            from uuid import UUID
            payload = {
                "source_type": src_type,
                "source_id": str(UUID(src_id.strip())),
                "target_type": tgt_type,
                "target_id": str(UUID(tgt_id.strip())),
            }
        except Exception:
            st.error("Invalid UUID format")
            return
        
        if rel_filter:
            payload["relationships"] = rel_filter
        
        with st.spinner("Finding path..."):
            try:
                result = _api_post("/api/graph/path", payload)
                
                if not result.get("found"):
                    st.info("No path found (or timed out)")
                else:
                    path_len = max(0, len(result.get("nodes", [])) - 1)
                    st.success(f"Path found! Length: {path_len} hops")
                    
                    # Show path nodes
                    with st.expander("Path Details"):
                        for i, node in enumerate(result.get("nodes", [])):
                            arrow = " → " if i < len(result.get("nodes", [])) - 1 else ""
                            st.write(f"{node.get('entity_type')}: {node.get('entity_id')}{arrow}")
            except Exception as e:
                st.error(f"Path finding failed: {e}")


def _render_compound_tab(layout: str, detect_clusters: bool, show_edge_labels: bool) -> None:
    """Compound-Target Networks tab."""
    st.subheader("Compound-Target Networks")
    st.info("🚧 Compound-Target network integration coming soon.")
    st.markdown("""
    **Features planned:**
    - Compound selector from database
    - Target interaction networks
    - Activity data overlays
    - Integration with existing compound-target component
    
    **For now**, use the existing [Compound-Target Network](/Compound%20Target%20Network) page.
    """)


def _render_pathways_tab() -> None:
    """Pathways tab - links to existing viewer."""
    st.subheader("Pathway Visualization")
    st.markdown("""
    **KEGG Pathway Maps** are available in the dedicated Pathway Map Viewer:
    - Browse and search KEGG pathways
    - Interactive pathway diagrams
    - Gene expression overlays
    - Data export capabilities
    """)
    
    if st.button("Open Pathway Map Viewer →"):
        st.switch_page("pages/pathway_map_viewer.py")


# Entry point for Streamlit multipage
def render() -> None:
    render_network_hub_page()


if __name__ == "__main__":
    render_network_hub_page()
