"""Unified Network Hub for Cytoscape.js visualizations."""

from __future__ import annotations

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
        _render_evidence_tab()
    
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


def _render_evidence_tab() -> None:
    """Evidence Graph tab - placeholder for migration."""
    st.subheader("Evidence Graph Explorer")
    st.info("🚧 Full Evidence Graph functionality will be migrated here in the next update.")
    st.markdown("""
    **Features coming soon:**
    - Entity traversal (k-hop neighborhoods)
    - Shortest path finding
    - Community detection visualization
    - Expression overlays on evidence networks
    
    **For now**, use the existing [Graph Explorer](/Graph%20Explorer) page.
    """)
    
    if st.button("Open Graph Explorer →"):
        st.switch_page("pages/graph_explorer.py")


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
