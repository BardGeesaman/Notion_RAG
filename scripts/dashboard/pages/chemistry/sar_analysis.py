"""SAR analysis tab."""

from __future__ import annotations

import math
import os
from typing import Dict, Any, Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import httpx


def ic50_to_pic50(ic50_nm: float) -> Optional[float]:
    """Convert IC50 (nM) to pIC50."""
    if ic50_nm is None or ic50_nm <= 0:
        return None
    return -math.log10(ic50_nm * 1e-9)


def _api_get(path: str) -> Dict[str, Any]:
    """Make GET request to API."""
    api_base = os.environ.get("API_URL", "http://localhost:8000")
    try:
        response = httpx.get(f"{api_base}{path}", timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def _api_post(path: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Make POST request to API."""
    api_base = os.environ.get("API_URL", "http://localhost:8000")
    try:
        response = httpx.post(f"{api_base}{path}", json=payload, timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def render_sar_tab(
    tab,
    *,
    db_session,
    compound_model,
    get_compound_properties,
    calculate_lipinski,
    get_activity_data,
    detect_activity_cliffs,
    find_common_core,
    decompose_rgroups,
    get_rgroup_statistics,
) -> None:
    """Render the SAR Analysis tab."""
    with tab:
        st.subheader("SAR Analysis")
        
        # Create sub-tabs for different SAR analyses
        sar_subtabs = st.tabs(["Property Analysis", "R-Group Table", "SAR Grid"])
        
        # Sub-tab 1: Property Analysis (existing content)
        with sar_subtabs[0]:
            _render_property_analysis(
                get_compound_properties, calculate_lipinski, 
                get_activity_data, detect_activity_cliffs
            )
        
        # Sub-tab 2: R-Group Table (existing R-group content)
        with sar_subtabs[1]:
            _render_rgroup_table(
                db_session, compound_model, find_common_core, 
                decompose_rgroups, get_rgroup_statistics
            )
        
        # Sub-tab 3: SAR Grid (new functionality)
        with sar_subtabs[2]:
            _render_sar_grid()


def _render_property_analysis(get_compound_properties, calculate_lipinski, get_activity_data, detect_activity_cliffs):
    """Render the property analysis sub-tab."""
    props_df = get_compound_properties()
    if props_df.empty:
        st.info("No compounds available for SAR analysis.")
        return

    lipinski_results = props_df.apply(lambda r: calculate_lipinski(r.get("smiles", "")), axis=1)
    props_df["lipinski_compliant"] = [res.get("passes_ro5", False) for res in lipinski_results]
    compliant = int(props_df["lipinski_compliant"].sum())
    st.metric("Compounds", len(props_df))
    st.metric("Lipinski Compliant", compliant)

    col_a, col_b = st.columns(2)
    with col_a:
        fig_mw = px.histogram(props_df, x="molecular_weight", nbins=40, title="Molecular Weight")
        st.plotly_chart(fig_mw, use_container_width=True)
    with col_b:
        fig_logp = px.histogram(props_df, x="logp", nbins=40, title="LogP")
        st.plotly_chart(fig_logp, use_container_width=True)

    act_df = get_activity_data()
    if not act_df.empty:
        merged = act_df.merge(
            props_df[["compound_id", "logp"]],
            on="compound_id",
            how="left",
        )
        merged_clean = merged.dropna(subset=["molecular_weight", "logp"])
        if not merged_clean.empty:
            scatter = px.scatter(
                merged_clean,
                x="logp",
                y="molecular_weight",
                title="LogP vs Molecular Weight (HTS)",
                hover_data=["compound_id", "smiles"],
            )
            st.plotly_chart(scatter, use_container_width=True)

        cliff_thresh = st.slider("Cliff similarity threshold", min_value=0.5, max_value=1.0, value=0.8, step=0.01)
        cliffs = detect_activity_cliffs(similarity_threshold=cliff_thresh)
        if cliffs:
            st.markdown("**Activity Cliffs (top 50 by activity difference)**")
            st.dataframe(cliffs[:50], hide_index=True, use_container_width=True)
        else:
            st.info("No activity cliffs detected at the selected threshold.")
    else:
        st.info("No HTS activity data available for cliff detection.")

    st.markdown("### Activity Cliffs Network")
    cliffs = detect_activity_cliffs()
    if cliffs:
        nodes = {}
        for c in cliffs:
            nodes[c["compound_1"]] = nodes.get(c["compound_1"], {"id": c["compound_1"]})
            nodes[c["compound_2"]] = nodes.get(c["compound_2"], {"id": c["compound_2"]})

        ids = list(nodes.keys())
        n = len(ids)
        for idx, cid in enumerate(ids):
            angle = 2 * math.pi * idx / max(n, 1)
            nodes[cid]["x"] = math.cos(angle)
            nodes[cid]["y"] = math.sin(angle)

        edge_x, edge_y, edge_text = [], [], []
        for c in cliffs:
            a, b = c["compound_1"], c["compound_2"]
            x0, y0 = nodes[a]["x"], nodes[a]["y"]
            x1, y1 = nodes[b]["x"], nodes[b]["y"]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
            edge_text.append(
                f"{a} ↔ {b}<br>similarity={c['similarity']:.2f}<br>fold_change={c['fold_change']:.1f}"
            )

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=1, color="lightgray"),
            hoverinfo="text",
            text=edge_text,
            mode="lines",
        )

        node_x = [nodes[cid]["x"] for cid in ids]
        node_y = [nodes[cid]["y"] for cid in ids]
        node_text = ids

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=node_text,
            textposition="bottom center",
            marker=dict(size=12, color="tomato"),
            hoverinfo="text",
        )

        fig = go.Figure(data=[edge_trace, node_trace])
        fig.update_layout(
            showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, visible=False),
            margin=dict(l=10, r=10, t=10, b=10),
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("No activity cliffs detected.")


def _render_rgroup_table(db_session, compound_model, find_common_core, decompose_rgroups, get_rgroup_statistics):
    """Render the R-group table sub-tab."""
    st.markdown("### R-Group Decomposition")
    st.markdown("Find common core and decompose compounds into R-groups for SAR analysis.")

    with db_session() as db:
        compounds = db.query(compound_model).limit(100).all()
        if compounds:
            compound_options = {f"{c.compound_id} - {c.smiles[:30]}": c.smiles for c in compounds}
            selected_compounds = st.multiselect(
                "Select compounds for R-group analysis",
                list(compound_options.keys()),
                help="Select 2 or more compounds to find common core",
            )

            if len(selected_compounds) >= 2:
                selected_smiles = [compound_options[c] for c in selected_compounds]

                col1, col2 = st.columns(2)
                with col1:
                    if st.button("Find Common Core", type="primary"):
                        with st.spinner("Finding maximum common substructure..."):
                            core_smarts = find_common_core(selected_smiles)
                            if core_smarts:
                                st.session_state["rgroup_core"] = core_smarts
                                st.success(f"Common core found: {core_smarts}")
                            else:
                                st.error("Could not find common core. Try different compounds.")

                if st.session_state.get("rgroup_core"):
                    core_smarts = st.session_state["rgroup_core"]
                    st.markdown(f"**Core SMARTS:** `{core_smarts}`")

                    with col2:
                        if st.button("Decompose R-Groups", type="primary"):
                            with st.spinner("Decomposing compounds..."):
                                decomposition = decompose_rgroups(selected_smiles, core_smarts)
                                if decomposition:
                                    st.session_state["rgroup_decomp"] = decomposition
                                    st.success(f"Decomposed {len(decomposition)} compounds")
                                else:
                                    st.error("R-group decomposition failed")

                    if st.session_state.get("rgroup_decomp"):
                        decomp = st.session_state["rgroup_decomp"]

                        st.markdown("**Decomposition Results:**")
                        decomp_df = pd.DataFrame(decomp)
                        st.dataframe(decomp_df, hide_index=True, use_container_width=True)

                        stats = get_rgroup_statistics(decomp)
                        if stats:
                            st.markdown("**R-Group Statistics:**")
                            for position, counts in stats.items():
                                st.markdown(f"**{position}:**")
                                stats_df = pd.DataFrame(
                                    [
                                        {"R-Group": smiles, "Frequency": count}
                                        for smiles, count in sorted(counts.items(), key=lambda x: x[1], reverse=True)
                                    ]
                                )
                                st.dataframe(stats_df, hide_index=True, use_container_width=True)
            else:
                st.info("Select at least 2 compounds to perform R-group analysis.")
        else:
            st.info("No compounds available for R-group analysis.")


def _render_sar_grid():
    """Render the SAR Grid sub-tab."""
    st.markdown("### SAR Grid Analysis")
    st.markdown("Generate structure-activity relationship heatmaps from R-group decomposition.")
    
    # Step 1: Target Selection
    st.markdown("#### 1. Target Selection")
    targets_data = _api_get("/api/v1/sar/targets")
    if not targets_data:
        st.error("Could not load targets from API")
        return
    
    target_options = [f"{t['target']} ({t['compound_count']} compounds)" for t in targets_data]
    if not target_options:
        st.info("No targets with compounds available")
        return
    
    selected_target = st.selectbox("Select target", target_options)
    if not selected_target:
        return
    
    # Extract target name
    target_name = selected_target.split(" (")[0]
    
    # Step 2: Fetch compounds for target
    compounds_data = _api_get(f"/api/v1/sar/targets/{target_name}/compounds")
    if not compounds_data:
        st.info("No compounds with activity data available for this target")
        return
    
    st.markdown(f"**Found {len(compounds_data)} compounds for target {target_name}**")
    
    # Step 3: Core Detection
    st.markdown("#### 2. Core Detection")
    compound_ids = [c["compound_id"] for c in compounds_data]
    
    col1, col2 = st.columns(2)
    with col1:
        auto_detect = st.button("Auto-detect Core", type="primary")
        
    with col2:
        manual_core = st.text_input("Manual Core SMARTS (optional)", placeholder="c1ccccc1")
    
    core_smiles = None
    if auto_detect:
        with st.spinner("Detecting common core..."):
            # Use a subset for core detection (performance)
            subset_ids = compound_ids[:20]  # Limit to first 20 compounds
            grid_request = {
                "compound_ids": subset_ids,
                "core_smarts": manual_core if manual_core else None,
                "x_axis": "R1",
                "y_axis": "R2"
            }
            
            grid_response = _api_post("/api/v1/sar/grid", grid_request)
            if grid_response and "core_smiles" in grid_response:
                core_smiles = grid_response["core_smiles"]
                st.session_state["sar_core"] = core_smiles
                st.success(f"Core detected: {core_smiles}")
            else:
                st.warning("Could not detect common core. Try manual SMARTS input.")
    
    if manual_core:
        core_smiles = manual_core
        st.session_state["sar_core"] = core_smiles
        st.info(f"Using manual core: {core_smiles}")
    
    if st.session_state.get("sar_core"):
        core_smiles = st.session_state["sar_core"]
        
        # Step 4: R-Group Axis Selection & Grid Generation
        st.markdown("#### 3. SAR Grid Configuration")
        
        col1, col2 = st.columns(2)
        with col1:
            x_axis = st.selectbox("X-Axis (columns)", ["R1", "R2", "R3", "R4"], index=0)
        with col2:
            y_axis = st.selectbox("Y-Axis (rows)", ["R1", "R2", "R3", "R4"], index=1)
        
        if x_axis == y_axis:
            st.warning("X-axis and Y-axis must be different")
            return
        
        # Step 5: Generate SAR Grid
        if st.button("Generate SAR Grid", type="primary"):
            with st.spinner("Generating SAR matrix..."):
                grid_request = {
                    "compound_ids": compound_ids,
                    "core_smarts": core_smiles,
                    "x_axis": x_axis,
                    "y_axis": y_axis
                }
                
                grid_response = _api_post("/api/v1/sar/grid", grid_request)
                
                if grid_response and "matrix" in grid_response:
                    st.session_state["sar_grid_data"] = grid_response
                    st.success(f"Generated SAR grid with {grid_response['total_compounds']} compounds")
                else:
                    st.error("Failed to generate SAR grid")
        
        # Step 6: Display Results
        if st.session_state.get("sar_grid_data"):
            grid_data = st.session_state["sar_grid_data"]
            
            # Convert to DataFrame for easier handling
            matrix = grid_data["matrix"]
            row_labels = grid_data["row_labels"]
            col_labels = grid_data["col_labels"]
            
            if not matrix or not row_labels or not col_labels:
                st.warning("Empty SAR matrix. Try different R-group positions.")
                return
            
            # Create DataFrame
            sar_df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
            
            # Convert IC50 to pIC50 for better visualization
            sar_pic50_df = sar_df.applymap(lambda x: ic50_to_pic50(x) if x is not None else None)
            
            # Step 7: Plotly Heatmap
            st.markdown("#### 4. SAR Heatmap")
            
            fig = px.imshow(
                sar_pic50_df.values,
                labels=dict(x=x_axis, y=y_axis, color="pIC50"),
                x=col_labels,
                y=row_labels,
                color_continuous_scale="RdYlGn",  # Red=low pIC50, Green=high pIC50
                aspect="auto",
                title=f"SAR Grid: {y_axis} vs {x_axis}"
            )
            
            fig.update_layout(
                xaxis_title=f"{x_axis} Substitutions",
                yaxis_title=f"{y_axis} Substitutions",
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Step 8: Statistics Summary
            st.markdown("#### 5. SAR Statistics")
            
            col1, col2, col3, col4 = st.columns(4)
            
            # Calculate statistics
            valid_values = [val for row in sar_pic50_df.values for val in row if val is not None]
            
            if valid_values:
                with col1:
                    st.metric(f"{x_axis} Variants", len(col_labels))
                with col2:
                    st.metric(f"{y_axis} Variants", len(row_labels))
                with col3:
                    best_pic50 = max(valid_values)
                    st.metric("Best pIC50", f"{best_pic50:.2f}")
                with col4:
                    flatness = pd.Series(valid_values).std() / pd.Series(valid_values).mean() if len(valid_values) > 1 else 0
                    st.metric("SAR Flatness", f"{flatness:.3f}")
                
                # Find best combination
                best_idx = sar_pic50_df.stack().idxmax()
                if best_idx:
                    st.info(f"Best combination: {y_axis}={best_idx[0]}, {x_axis}={best_idx[1]} (pIC50={best_pic50:.2f})")
            
            # Step 9: Download Button
            st.markdown("#### 6. Export Data")
            
            # Prepare CSV data
            csv_data = sar_df.copy()
            csv_data.index.name = y_axis
            csv_data.columns.name = x_axis
            
            csv_string = csv_data.to_csv()
            
            st.download_button(
                label="Download SAR Matrix (CSV)",
                data=csv_string,
                file_name=f"sar_matrix_{target_name}_{x_axis}_{y_axis}.csv",
                mime="text/csv"
            )

