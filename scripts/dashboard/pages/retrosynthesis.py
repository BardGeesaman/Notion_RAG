"""Retrosynthesis Advisor Dashboard - AI-powered synthesis route planning."""

from __future__ import annotations

import os
from typing import Dict, Optional

import httpx
import streamlit as st

from scripts.dashboard.auth import require_auth
from scripts.dashboard.components.cytoscape import render_cytoscape


def _api_base() -> str:
    return os.environ.get("API_URL", "http://localhost:8000").rstrip("/")


def _api_post(path: str, payload: Dict) -> Dict:
    try:
        with httpx.Client(timeout=60) as client:
            r = client.post(f"{_api_base()}{path}", json=payload)
            r.raise_for_status()
            return r.json()
    except Exception as e:
        st.error(f"API error: {e}")
        return {}


def _api_get(path: str, params: Optional[Dict] = None) -> Dict:
    try:
        with httpx.Client(timeout=60) as client:
            r = client.get(f"{_api_base()}{path}", params=params)
            r.raise_for_status()
            return r.json()
    except Exception as e:
        st.error(f"API error: {e}")
        return {}


def render_retrosynthesis_page() -> None:
    """Render the Retrosynthesis Advisor page."""
    require_auth()
    
    st.header("🧪 Retrosynthesis Advisor")
    st.caption("AI-powered synthesis route planning for target molecules.")
    
    # Session state
    if "retro_analysis" not in st.session_state:
        st.session_state.retro_analysis = None
    
    # Tabs
    tab1, tab2, tab3 = st.tabs(["🎯 Analyze", "🔀 Routes", "🧱 Building Blocks"])
    
    with tab1:
        render_analyze_tab()
    
    with tab2:
        render_routes_tab()
    
    with tab3:
        render_building_blocks_tab()


def render_analyze_tab() -> None:
    """Render the analysis tab with SMILES input and tree visualization."""
    st.subheader("Target Analysis")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("### Input")
        
        smiles_input = st.text_area(
            "Target SMILES",
            placeholder="Enter target molecule SMILES\nExample: CC(=O)Nc1ccc(C(=O)O)cc1",
            height=100,
            key="target_smiles"
        )
        
        max_depth = st.slider(
            "Max Synthesis Depth",
            min_value=1,
            max_value=10,
            value=5,
            help="Maximum number of synthetic steps"
        )
        
        if st.button("🔍 Analyze Target", type="primary", disabled=not smiles_input.strip()):
            with st.spinner("Analyzing synthesis routes..."):
                result = _api_post("/api/v1/retrosynthesis/analyze", {
                    "smiles": smiles_input.strip(),
                    "max_depth": max_depth
                })
                
                if result:
                    st.session_state.retro_analysis = result
                    st.success(f"Found {result.get('tree', {}).get('num_alternatives', 0)} synthesis routes!")
    
    with col2:
        st.markdown("### Synthesis Tree")
        
        analysis = st.session_state.retro_analysis
        if analysis and analysis.get("tree"):
            tree = analysis["tree"]
            
            # Convert to Cytoscape format
            nodes, edges = synthesis_tree_to_cytoscape(tree)
            
            if nodes and edges:
                render_cytoscape(
                    nodes=nodes,
                    edges=edges,
                    height=500,
                    layout="breadthfirst",
                    show_edge_labels=True,
                    show_legend=False
                )
            else:
                st.info("Tree visualization will appear here after analysis")
        else:
            st.info("Enter a target SMILES and click Analyze to see synthesis routes")


def render_routes_tab() -> None:
    """Render route comparison tab."""
    st.subheader("Route Comparison")
    
    analysis = st.session_state.retro_analysis
    if not analysis or not analysis.get("tree"):
        st.info("🎯 Run analysis first to see route comparisons")
        return
    
    tree = analysis["tree"]
    routes = tree.get("routes", [])
    
    if not routes:
        st.warning("No routes found in analysis")
        return
    
    # Route selector
    route_options = {}
    for i, route in enumerate(routes):
        # Calculate score for display
        score = _calculate_display_score(route)
        route_options[f"Route {i+1} ({score:.1f} pts, {route['total_steps']} steps)"] = i
    
    selected_route_name = st.selectbox(
        "Select Route",
        options=list(route_options.keys()),
        key="selected_route"
    )
    
    if selected_route_name:
        route_idx = route_options[selected_route_name]
        route = routes[route_idx]
        
        # Score breakdown
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("### Route Details")
            st.metric("Total Steps", route["total_steps"])
            st.metric("Confidence", f"{route['confidence']:.2f}")
            
            # Steps table
            st.markdown("### Synthesis Steps")
            for i, step in enumerate(route["steps"]):
                with st.expander(f"Step {i+1}: {step['reaction_type']}"):
                    st.write(f"**Reactants:** {', '.join(step['reactants'])}")
                    st.write(f"**Product:** {step['product']}")
                    st.write(f"**Conditions:** {step['conditions']}")
                    st.write(f"**Confidence:** {step['confidence']:.2f}")
        
        with col2:
            st.markdown("### Route Scoring")
            
            # Get detailed score
            score_result = _api_post("/api/v1/retrosynthesis/score", route)
            if score_result:
                st.metric("Total Score", f"{score_result['total_score']:.1f}/100")
                st.metric("Complexity Score", f"{score_result['complexity_score']:.1f}")
                st.metric("Availability Score", f"{score_result['availability_score']:.1f}")
                st.metric("Estimated Cost", f"${score_result['estimated_cost']:.0f}")


def render_building_blocks_tab() -> None:
    """Render building block availability checker."""
    st.subheader("Building Block Availability")
    
    smiles_input = st.text_area(
        "Building Block SMILES",
        placeholder="Enter SMILES strings (one per line)\nExample:\nCC(C)OC(=O)c1ccc(N)cc1\nCC(=O)Cl",
        height=150,
        key="bb_smiles"
    )
    
    if st.button("🔍 Check Availability", disabled=not smiles_input.strip()):
        smiles_list = [s.strip() for s in smiles_input.strip().split("\n") if s.strip()]
        
        if smiles_list:
            with st.spinner("Checking building block availability..."):
                # Note: API expects query params, but using POST for multiple SMILES
                results = []
                for smiles in smiles_list:
                    result = _api_get("/api/v1/retrosynthesis/building-blocks", {"smiles": [smiles]})
                    if result:
                        results.extend(result)
                
                if results:
                    # Display results table
                    st.markdown("### Availability Results")
                    
                    for result in results:
                        col1, col2, col3 = st.columns([2, 1, 2])
                        
                        with col1:
                            st.code(result["smiles"])
                        
                        with col2:
                            if result["available"]:
                                st.success("✅ Available")
                            else:
                                st.error("❌ Not Available")
                        
                        with col3:
                            if result.get("vendors"):
                                vendor_text = ", ".join([v["name"] for v in result["vendors"]])
                                st.write(vendor_text)
                            else:
                                st.write("No vendors")


def synthesis_tree_to_cytoscape(tree: dict) -> tuple[list, list]:
    """Convert synthesis tree to Cytoscape nodes and edges."""
    nodes = []
    edges = []
    
    target_smiles = tree["target"]
    
    # Target node (green)
    target_id = "target"
    nodes.append({
        "data": {
            "id": target_id,
            "label": target_smiles[:20] + ("..." if len(target_smiles) > 20 else ""),
            "smiles": target_smiles,
            "node_type": "target",
            "color": "#28a745"  # Green
        }
    })
    
    # Process routes to create tree structure
    routes = tree.get("routes", [])
    
    for route_idx, route in enumerate(routes):
        steps = route.get("steps", [])
        
        # Create nodes for each step
        for step_idx, step in enumerate(steps):
            # Reactant nodes
            for reactant_idx, reactant in enumerate(step["reactants"]):
                reactant_id = f"route{route_idx}_step{step_idx}_reactant{reactant_idx}"
                
                nodes.append({
                    "data": {
                        "id": reactant_id,
                        "label": reactant[:15] + ("..." if len(reactant) > 15 else ""),
                        "smiles": reactant,
                        "node_type": "building_block",
                        "color": "#6c757d"  # Gray
                    }
                })
                
                # Edge from reactant to product
                product_id = f"route{route_idx}_step{step_idx}_product"
                if step_idx == len(steps) - 1:
                    # Last step connects to target
                    product_id = target_id
                
                edges.append({
                    "data": {
                        "id": f"{reactant_id}_to_{product_id}",
                        "source": reactant_id,
                        "target": product_id,
                        "label": step["reaction_type"],
                        "confidence": step["confidence"]
                    }
                })
            
            # Intermediate product nodes (except final target)
            if step_idx < len(steps) - 1:
                product_id = f"route{route_idx}_step{step_idx}_product"
                nodes.append({
                    "data": {
                        "id": product_id,
                        "label": step["product"][:15] + ("..." if len(step["product"]) > 15 else ""),
                        "smiles": step["product"],
                        "node_type": "intermediate",
                        "color": "#007bff"  # Blue
                    }
                })
    
    return nodes, edges


def _calculate_display_score(route: dict) -> float:
    """Calculate a simple display score for route selection."""
    # Simple scoring for display
    confidence = route.get("confidence", 0.5)
    step_penalty = route.get("total_steps", 5) * 5  # Fewer steps = higher score
    return max(0, confidence * 100 - step_penalty)


# Page registration
def render() -> None:
    render_retrosynthesis_page()


if __name__ == "__main__":
    render_retrosynthesis_page()
