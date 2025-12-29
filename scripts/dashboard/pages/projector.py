"""High-Dimensional Projector page."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_projector_page() -> None:
    """Render the High-Dimensional Projector page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📊 High-Dimensional Projector")
    st.caption("Visualize high-dimensional data using UMAP, t-SNE, or PCA.")
    
    # Sidebar controls
    with st.sidebar:
        st.subheader("Projection Settings")
        
        # Dataset selector
        if st.button("Load Datasets"):
            try:
                response = _api_get("/api/v1/projector/datasets")
                st.session_state["available_datasets"] = response.get("datasets", [])
            except Exception as e:
                st.error(f"Failed to load datasets: {e}")
        
        datasets = st.session_state.get("available_datasets", [])
        
        if datasets:
            dataset_options = {f"{ds['name']} ({ds['omics_type']})": ds['id'] for ds in datasets}
            selected_name = st.selectbox("Select Dataset", list(dataset_options.keys()))
            selected_dataset_id = dataset_options[selected_name]
        else:
            st.info("Click 'Load Datasets' to view available data")
            selected_dataset_id = None
        
        # Algorithm selector
        algorithm = st.selectbox(
            "Algorithm",
            options=["UMAP", "t-SNE", "PCA"],
            help="Dimensionality reduction algorithm",
        )
        
        # Color-by selector
        color_by = st.selectbox(
            "Color By",
            options=["None", "Index"],
            help="Color points by metadata",
        )
        
        # Dimensionality
        n_components = st.radio("Dimensions", options=[2, 3], horizontal=True, index=0)
        
        # Algorithm-specific parameters
        if algorithm == "UMAP":
            n_neighbors = st.slider("n_neighbors", min_value=2, max_value=100, value=15, step=1)
            min_dist = st.slider("min_dist", min_value=0.0, max_value=1.0, value=0.1, step=0.05)
        elif algorithm == "t-SNE":
            perplexity = st.slider("Perplexity", min_value=5, max_value=50, value=30, step=5)
        
        random_state = st.number_input("Random Seed", min_value=0, value=42, step=1)
        
        compute_button = st.button("🚀 Compute Projection", type="primary", disabled=not selected_dataset_id)
        
        # Store color setting
        st.session_state["color_by_setting"] = color_by
    
    # Main area
    if compute_button and selected_dataset_id:
        with st.spinner(f"Computing {algorithm} projection..."):
            try:
                payload = {
                    "dataset_id": selected_dataset_id,
                    "algorithm": algorithm.lower(),
                    "n_components": n_components,
                    "random_state": random_state,
                }
                
                if algorithm == "UMAP":
                    payload["n_neighbors"] = n_neighbors
                    payload["min_dist"] = min_dist
                elif algorithm == "t-SNE":
                    payload["perplexity"] = perplexity
                
                result = _api_post("/api/v1/projector/compute", payload, timeout=180)
                st.session_state["projection_result"] = result
            except Exception as e:
                st.error(f"Projection failed: {e}")
                return
    
    # Display projection
    result = st.session_state.get("projection_result")
    
    if result:
        coords = result.get("coordinates", [])
        n_samples = result.get("n_samples", 0)
        algorithm_used = result.get("algorithm_used", "unknown")
        cached = result.get("cached", False)
        
        st.success(f"Projected {n_samples} samples using {algorithm_used.upper()} {'(cached)' if cached else ''}")
        
        # Create Plotly visualization
        color_by_setting = st.session_state.get("color_by_setting", "None")
        
        if n_components == 2:
            import plotly.express as px
            import pandas as pd
            
            df = pd.DataFrame(coords, columns=["X", "Y"])
            df["Sample"] = [f"sample_{i}" for i in range(len(coords))]
            
            # Add color column if requested
            color_col = None
            if color_by_setting == "Index":
                df["Color"] = df.index
                color_col = "Color"
            
            fig = px.scatter(
                df,
                x="X",
                y="Y",
                color=color_col,
                hover_data=["Sample"],
                title=f"{algorithm_used.upper()} 2D Projection",
            )
            st.plotly_chart(fig, use_container_width=True)
        
        elif n_components == 3:
            import plotly.express as px
            import pandas as pd
            
            df = pd.DataFrame(coords, columns=["X", "Y", "Z"])
            df["Sample"] = [f"sample_{i}" for i in range(len(coords))]
            
            # Add color column if requested
            color_col = None
            if color_by_setting == "Index":
                df["Color"] = df.index
                color_col = "Color"
            
            fig = px.scatter_3d(
                df,
                x="X",
                y="Y",
                z="Z",
                color=color_col,
                hover_data=["Sample"],
                title=f"{algorithm_used.upper()} 3D Projection",
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # Export button
        if st.button("📥 Export Coordinates (CSV)"):
            try:
                export_result = _api_post("/api/v1/projector/export", {"coordinates": coords})
                csv_data = export_result.get("csv_data", "")
                filename = export_result.get("filename", "projection.csv")
                
                st.download_button(
                    "Download CSV",
                    data=csv_data,
                    file_name=filename,
                    mime="text/csv",
                )
            except Exception as e:
                st.error(f"Export failed: {e}")
    
    else:
        st.info("Select a dataset and algorithm from the sidebar, then click 'Compute Projection'.")


if __name__ == "__main__":
    render_projector_page()

