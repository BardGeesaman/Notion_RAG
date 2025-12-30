"""Flow Cytometry dashboard page."""

from __future__ import annotations

import os
import tempfile
from typing import Any, Dict, List, Optional
import json

import httpx
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_flow_cytometry_page() -> None:
    """Render the Flow Cytometry analysis dashboard."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🔬 Flow Cytometry Analysis")
    st.caption("Upload FCS files, create gates, and analyze cell populations.")

    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📤 Upload", "📊 Scatter Plots", "🎯 Gating", "📈 Statistics"])

    with tab1:
        render_upload_tab()
    
    with tab2:
        render_scatter_plots_tab()
    
    with tab3:
        render_gating_tab()
    
    with tab4:
        render_statistics_tab()


def render_upload_tab() -> None:
    """Render the Upload tab for FCS file upload."""
    st.subheader("Upload FCS Files")
    
    # File uploader
    uploaded_file = st.file_uploader(
        "Choose an FCS file",
        type=['fcs'],
        help="Upload Flow Cytometry Standard (.fcs) files for analysis"
    )
    
    if uploaded_file is not None:
        st.success(f"File selected: {uploaded_file.name}")
        
        # Show file details
        col1, col2 = st.columns(2)
        with col1:
            st.metric("File Size", f"{len(uploaded_file.getvalue()) / 1024:.1f} KB")
        with col2:
            st.metric("File Type", uploaded_file.type or "application/octet-stream")
        
        # Upload button
        if st.button("Upload and Process", type="primary", use_container_width=True):
            with st.spinner("Uploading and processing FCS file..."):
                try:
                    # Create temporary file
                    with tempfile.NamedTemporaryFile(suffix=".fcs", delete=False) as tmp_file:
                        tmp_file.write(uploaded_file.getvalue())
                        tmp_file.flush()
                        
                        # Upload via API
                        with httpx.Client(timeout=120) as client:
                            with open(tmp_file.name, "rb") as f:
                                files = {"file": (uploaded_file.name, f, "application/octet-stream")}
                                response = client.post(f"{API_BASE}/api/v1/flow-cytometry/upload", files=files)
                        
                        if response.status_code == 201:
                            data = response.json()
                            st.success("✅ File uploaded successfully!")
                            
                            # Show processing status
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Dataset ID", data["flow_dataset_id"][:8] + "...")
                            with col2:
                                st.metric("Processing Status", data["processing_status"])
                            with col3:
                                st.metric("File Size", f"{data['file_size_bytes'] / 1024:.1f} KB")
                            
                            # Store dataset ID in session state for other tabs
                            st.session_state.current_dataset_id = data["flow_dataset_id"]
                            
                            st.info("💡 Switch to the 'Scatter Plots' tab to visualize your data!")
                            
                        else:
                            st.error(f"Upload failed: {response.text}")
                            
                except Exception as e:
                    st.error(f"Error uploading file: {str(e)}")
    
    # Show existing datasets
    st.subheader("Existing Datasets")
    
    try:
        with st.spinner("Loading datasets..."):
            with httpx.Client(timeout=30) as client:
                response = client.get(f"{API_BASE}/api/v1/flow-cytometry/datasets", params={"limit": 50})
        
        if response.status_code == 200:
            data = response.json()
            datasets = data.get("items", [])
            
            if datasets:
                # Create a dataframe for display
                df_datasets = pd.DataFrame([
                    {
                        "ID": d["id"][:8] + "...",
                        "Status": d["processing_status"],
                        "Events": f"{d.get('n_events', 0):,}" if d.get('n_events') else "N/A",
                        "Parameters": d.get('n_parameters', 0),
                        "Cytometer": d.get('cytometer_model', 'Unknown'),
                        "Sample": d.get('sample_id', 'N/A'),
                    }
                    for d in datasets
                ])
                
                # Dataset selector
                selected_idx = st.selectbox(
                    "Select a dataset to analyze",
                    range(len(datasets)),
                    format_func=lambda i: f"{datasets[i]['id'][:8]}... - {datasets[i]['processing_status']} - {datasets[i].get('n_events', 0):,} events",
                    help="Choose a dataset to work with in other tabs"
                )
                
                if selected_idx is not None:
                    selected_dataset = datasets[selected_idx]
                    st.session_state.current_dataset_id = selected_dataset["id"]
                    st.session_state.current_dataset = selected_dataset
                
                # Show dataset table
                st.dataframe(df_datasets, use_container_width=True)
                
            else:
                st.info("No datasets found. Upload an FCS file to get started!")
        else:
            st.error(f"Failed to load datasets: {response.text}")
            
    except Exception as e:
        st.error(f"Error loading datasets: {str(e)}")


def render_scatter_plots_tab() -> None:
    """Render the Scatter Plots tab for data visualization."""
    st.subheader("2D Scatter Plots")
    
    # Check if dataset is selected
    if "current_dataset_id" not in st.session_state:
        st.warning("Please select a dataset in the Upload tab first.")
        return
    
    dataset_id = st.session_state.current_dataset_id
    
    # Load parameters
    try:
        with st.spinner("Loading parameters..."):
            with httpx.Client(timeout=30) as client:
                response = client.get(f"{API_BASE}/api/v1/flow-cytometry/datasets/{dataset_id}/parameters")
        
        if response.status_code == 200:
            parameters = response.json()
            param_names = [p["parameter_name"] for p in parameters]
            
            if len(param_names) < 2:
                st.error("Dataset needs at least 2 parameters for scatter plots.")
                return
            
            # Parameter selection
            col1, col2 = st.columns(2)
            with col1:
                x_param = st.selectbox("X Parameter", param_names, index=0)
            with col2:
                y_param = st.selectbox("Y Parameter", param_names, index=1 if len(param_names) > 1 else 0)
            
            # Transformation options
            col1, col2, col3 = st.columns(3)
            with col1:
                x_transform = st.selectbox("X Transform", ["linear", "log", "logicle", "arcsinh"], key="x_transform")
            with col2:
                y_transform = st.selectbox("Y Transform", ["linear", "log", "logicle", "arcsinh"], key="y_transform")
            with col3:
                density_mode = st.checkbox("Density coloring", value=True, help="Use hexbin density coloring for better visualization of dense regions")
            
            # Load events data
            if st.button("Generate Plot", type="primary", use_container_width=True):
                with st.spinner("Loading event data..."):
                    try:
                        with httpx.Client(timeout=60) as client:
                            response = client.get(
                                f"{API_BASE}/api/v1/flow-cytometry/datasets/{dataset_id}/events",
                                params={"subsample": True, "limit": 50000}
                            )
                        
                        if response.status_code == 200:
                            data = response.json()
                            events = data["events"]
                            
                            if events:
                                # Convert to DataFrame
                                df = pd.DataFrame(events)
                                
                                # Apply transformations (simplified - normally would use proper flow transforms)
                                import numpy as np
                                if x_transform == "log":
                                    df[x_param] = df[x_param].apply(lambda x: max(1, x)).apply(lambda x: np.log10(x))
                                elif x_transform == "arcsinh":
                                    df[x_param] = df[x_param].apply(lambda x: np.arcsinh(x / 150))
                                
                                if y_transform == "log":
                                    df[y_param] = df[y_param].apply(lambda x: max(1, x)).apply(lambda x: np.log10(x))
                                elif y_transform == "arcsinh":
                                    df[y_param] = df[y_param].apply(lambda x: np.arcsinh(x / 150))
                                
                                # Create plot
                                if density_mode:
                                    fig = px.density_heatmap(
                                        df, x=x_param, y=y_param,
                                        title=f"{y_param} vs {x_param}",
                                        nbinsx=50, nbinsy=50,
                                        color_continuous_scale="Viridis"
                                    )
                                else:
                                    # Sample down for scatter plot to avoid overplotting
                                    if len(df) > 10000:
                                        df_sample = df.sample(n=10000)
                                    else:
                                        df_sample = df
                                    
                                    fig = px.scatter(
                                        df_sample, x=x_param, y=y_param,
                                        title=f"{y_param} vs {x_param}",
                                        opacity=0.6
                                    )
                                
                                fig.update_layout(
                                    width=700,
                                    height=500,
                                    xaxis_title=f"{x_param} ({x_transform})",
                                    yaxis_title=f"{y_param} ({y_transform})"
                                )
                                
                                st.plotly_chart(fig, use_container_width=True)
                                
                                # Show data summary
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Total Events", f"{data['total_events']:,}")
                                with col2:
                                    st.metric("Displayed Events", f"{len(events):,}")
                                with col3:
                                    st.metric("Parameters", len(df.columns))
                                
                                # Store plot data for gating
                                st.session_state.plot_data = df
                                st.session_state.current_x_param = x_param
                                st.session_state.current_y_param = y_param
                                
                            else:
                                st.warning("No events found in dataset.")
                        else:
                            st.error(f"Failed to load events: {response.text}")
                            
                    except Exception as e:
                        st.error(f"Error loading events: {str(e)}")
        else:
            st.error(f"Failed to load parameters: {response.text}")
            
    except Exception as e:
        st.error(f"Error loading parameters: {str(e)}")


def render_gating_tab() -> None:
    """Render the Gating tab for interactive gate creation."""
    st.subheader("Interactive Gating")
    
    # Check if dataset is selected
    if "current_dataset_id" not in st.session_state:
        st.warning("Please select a dataset in the Upload tab first.")
        return
    
    dataset_id = st.session_state.current_dataset_id
    
    # Gate creation form
    st.subheader("Create New Gate")
    
    gate_type = st.selectbox("Gate Type", ["polygon", "rectangle", "quadrant"], help="Select the type of gate to create")
    gate_name = st.text_input("Gate Name", placeholder="Enter a name for this gate")
    
    if gate_type == "polygon":
        st.info("💡 Polygon gates: Define vertices as coordinate pairs")
        vertices_text = st.text_area(
            "Vertices (JSON format)",
            placeholder='[[100, 100], [900, 100], [900, 900], [100, 900]]',
            help="Enter polygon vertices as a JSON array of [x, y] coordinates"
        )
        
        if gate_name and vertices_text:
            try:
                vertices = json.loads(vertices_text)
                gate_definition = {"vertices": vertices}
                
                # Get current parameters
                x_param = st.session_state.get("current_x_param", "FSC-A")
                y_param = st.session_state.get("current_y_param", "SSC-A")
                
                if st.button("Create Polygon Gate", type="primary"):
                    create_gate(dataset_id, gate_name, gate_type, x_param, y_param, gate_definition)
                    
            except json.JSONDecodeError:
                st.error("Invalid JSON format for vertices")
    
    elif gate_type == "rectangle":
        st.info("💡 Rectangle gates: Define min/max bounds for X and Y parameters")
        col1, col2 = st.columns(2)
        with col1:
            x_min = st.number_input("X Min", value=0.0)
            x_max = st.number_input("X Max", value=1000.0)
        with col2:
            y_min = st.number_input("Y Min", value=0.0)
            y_max = st.number_input("Y Max", value=1000.0)
        
        if gate_name:
            gate_definition = {"x_min": x_min, "x_max": x_max, "y_min": y_min, "y_max": y_max}
            x_param = st.session_state.get("current_x_param", "FSC-A")
            y_param = st.session_state.get("current_y_param", "SSC-A")
            
            if st.button("Create Rectangle Gate", type="primary"):
                create_gate(dataset_id, gate_name, gate_type, x_param, y_param, gate_definition)
    
    elif gate_type == "quadrant":
        st.info("💡 Quadrant gates: Define X and Y thresholds to create 4 regions")
        col1, col2 = st.columns(2)
        with col1:
            x_threshold = st.number_input("X Threshold", value=500.0)
        with col2:
            y_threshold = st.number_input("Y Threshold", value=500.0)
        
        if gate_name:
            gate_definition = {"x_thresh": x_threshold, "y_thresh": y_threshold}
            x_param = st.session_state.get("current_x_param", "FSC-A")
            y_param = st.session_state.get("current_y_param", "SSC-A")
            
            if st.button("Create Quadrant Gate", type="primary"):
                create_gate(dataset_id, gate_name, gate_type, x_param, y_param, gate_definition)
    
    # Show existing gates
    st.subheader("Existing Gates")
    
    try:
        with st.spinner("Loading gates..."):
            with httpx.Client(timeout=30) as client:
                response = client.get(f"{API_BASE}/api/v1/flow-cytometry/datasets/{dataset_id}/gates")
        
        if response.status_code == 200:
            gates = response.json()
            
            if gates:
                # Create gates dataframe
                df_gates = pd.DataFrame([
                    {
                        "Name": g["gate_name"],
                        "Type": g["gate_type"],
                        "X Parameter": g.get("x_parameter", "N/A"),
                        "Y Parameter": g.get("y_parameter", "N/A"),
                        "Active": "✅" if g["is_active"] else "❌",
                        "ID": g["id"][:8] + "..."
                    }
                    for g in gates
                ])
                
                st.dataframe(df_gates, use_container_width=True)
                
                # Gate management
                st.subheader("Gate Management")
                selected_gate_idx = st.selectbox(
                    "Select gate to modify",
                    range(len(gates)),
                    format_func=lambda i: gates[i]["gate_name"]
                )
                
                if selected_gate_idx is not None:
                    selected_gate = gates[selected_gate_idx]
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if st.button("Deactivate Gate", type="secondary"):
                            deactivate_gate(selected_gate["id"])
                    
                    with col2:
                        if st.button("Delete Gate", type="secondary"):
                            delete_gate(selected_gate["id"])
            else:
                st.info("No gates created yet. Create a gate above to get started!")
        else:
            st.error(f"Failed to load gates: {response.text}")
            
    except Exception as e:
        st.error(f"Error loading gates: {str(e)}")


def render_statistics_tab() -> None:
    """Render the Statistics tab for population analysis."""
    st.subheader("Population Statistics")
    
    # Check if dataset is selected
    if "current_dataset_id" not in st.session_state:
        st.warning("Please select a dataset in the Upload tab first.")
        return
    
    dataset_id = st.session_state.current_dataset_id
    
    # Load population statistics
    try:
        with st.spinner("Loading population statistics..."):
            with httpx.Client(timeout=30) as client:
                response = client.get(f"{API_BASE}/api/v1/flow-cytometry/datasets/{dataset_id}/populations")
        
        if response.status_code == 200:
            populations = response.json()
            
            if populations:
                # Create population statistics table
                df_populations = pd.DataFrame([
                    {
                        "Population": p["population_name"],
                        "Events": f"{p['n_events']:,}",
                        "% of Parent": f"{p.get('pct_of_parent', 0):.1f}%",
                        "% of Total": f"{p.get('pct_of_total', 0):.1f}%",
                        "Gate ID": p.get("gate_id", "N/A")[:8] + "..." if p.get("gate_id") else "N/A"
                    }
                    for p in populations
                ])
                
                st.dataframe(df_populations, use_container_width=True)
                
                # Summary metrics
                total_events = sum(p["n_events"] for p in populations)
                st.subheader("Summary")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Populations", len(populations))
                with col2:
                    st.metric("Total Gated Events", f"{total_events:,}")
                with col3:
                    st.metric("Average Population Size", f"{total_events // len(populations) if populations else 0:,}")
                
                # Population distribution pie chart
                if len(populations) > 1:
                    st.subheader("Population Distribution")
                    fig = px.pie(
                        values=[p["n_events"] for p in populations],
                        names=[p["population_name"] for p in populations],
                        title="Population Distribution by Event Count"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                # Export functionality
                st.subheader("Export Data")
                if st.button("Download Population Statistics", type="primary"):
                    csv_data = df_populations.to_csv(index=False)
                    st.download_button(
                        label="📥 Download CSV",
                        data=csv_data,
                        file_name=f"population_stats_{dataset_id[:8]}.csv",
                        mime="text/csv"
                    )
                
                # Channel statistics (if available)
                st.subheader("Channel Statistics")
                for pop in populations:
                    if pop.get("median_values") or pop.get("mean_values"):
                        with st.expander(f"📊 {pop['population_name']} Channel Stats"):
                            if pop.get("median_values"):
                                st.write("**Median Values:**")
                                median_df = pd.DataFrame(list(pop["median_values"].items()), columns=["Parameter", "Median"])
                                st.dataframe(median_df, use_container_width=True)
                            
                            if pop.get("mean_values"):
                                st.write("**Mean Values:**")
                                mean_df = pd.DataFrame(list(pop["mean_values"].items()), columns=["Parameter", "Mean"])
                                st.dataframe(mean_df, use_container_width=True)
                            
                            if pop.get("cv_values"):
                                st.write("**CV Values:**")
                                cv_df = pd.DataFrame(list(pop["cv_values"].items()), columns=["Parameter", "CV %"])
                                st.dataframe(cv_df, use_container_width=True)
                
            else:
                st.info("No population statistics available. Create gates in the Gating tab to generate populations!")
        else:
            st.error(f"Failed to load population statistics: {response.text}")
            
    except Exception as e:
        st.error(f"Error loading population statistics: {str(e)}")


def create_gate(dataset_id: str, gate_name: str, gate_type: str, x_parameter: str, y_parameter: str, gate_definition: dict) -> None:
    """Create a new gate via API."""
    try:
        gate_data = {
            "gate_name": gate_name,
            "gate_type": gate_type,
            "x_parameter": x_parameter,
            "y_parameter": y_parameter,
            "gate_definition": gate_definition
        }
        
        with httpx.Client(timeout=60) as client:
            response = client.post(
                f"{API_BASE}/api/v1/flow-cytometry/datasets/{dataset_id}/gates",
                json=gate_data
            )
        
        if response.status_code == 201:
            st.success(f"✅ Gate '{gate_name}' created successfully!")
            st.rerun()
        else:
            st.error(f"Failed to create gate: {response.text}")
            
    except Exception as e:
        st.error(f"Error creating gate: {str(e)}")


def deactivate_gate(gate_id: str) -> None:
    """Deactivate a gate via API."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.put(
                f"{API_BASE}/api/v1/flow-cytometry/gates/{gate_id}",
                json={"is_active": False}
            )
        
        if response.status_code == 200:
            st.success("✅ Gate deactivated successfully!")
            st.rerun()
        else:
            st.error(f"Failed to deactivate gate: {response.text}")
            
    except Exception as e:
        st.error(f"Error deactivating gate: {str(e)}")


def delete_gate(gate_id: str) -> None:
    """Delete a gate via API."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.delete(f"{API_BASE}/api/v1/flow-cytometry/gates/{gate_id}")
        
        if response.status_code == 204:
            st.success("✅ Gate deleted successfully!")
            st.rerun()
        else:
            st.error(f"Failed to delete gate: {response.text}")
            
    except Exception as e:
        st.error(f"Error deleting gate: {str(e)}")


if __name__ == "__main__":
    render_flow_cytometry_page()
