"""Active Learning Dashboard."""

import os
import uuid
from typing import Dict, Any, List, Optional

import pandas as pd
import streamlit as st
import requests
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime

# API Configuration
API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def get_auth_headers() -> Dict[str, str]:
    """Get authentication headers for API requests."""
    # In production, this would use actual session tokens
    return {"Authorization": "Bearer dummy-token"}


def api_get(path: str) -> Dict[str, Any]:
    """Make authenticated GET request to API."""
    try:
        response = requests.get(f"{API_BASE}{path}", headers=get_auth_headers(), timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def api_post(path: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Make authenticated POST request to API."""
    try:
        response = requests.post(f"{API_BASE}{path}", json=payload, headers=get_auth_headers(), timeout=30.0)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def get_uncertainty_color(uncertainty: Optional[float]) -> str:
    """Get color for uncertainty visualization."""
    if uncertainty is None:
        return "gray"
    elif uncertainty > 0.4:
        return "red"    # High uncertainty
    elif uncertainty > 0.2:
        return "orange" # Medium uncertainty
    else:
        return "green"  # Low uncertainty


def get_ml_models() -> List[Dict[str, Any]]:
    """Get list of available ML models."""
    # Mock for MVP - in production, query ML models API
    return [
        {"id": str(uuid.uuid4()), "name": "ADMET hERG Predictor", "model_type": "classification"},
        {"id": str(uuid.uuid4()), "name": "LogP Regression", "model_type": "regression"},
        {"id": str(uuid.uuid4()), "name": "Solubility Predictor", "model_type": "regression"},
    ]


def render_label_queue_tab():
    """Render Tab 1: Label Queue."""
    st.subheader("Label Queue")
    
    # Model selection
    models = get_ml_models()
    if not models:
        st.info("No ML models available for active learning.")
        return
    
    model_options = {f"{m['name']} ({m['model_type']})": m["id"] for m in models}
    selected_model_key = st.selectbox("Select Model", list(model_options.keys()))
    
    if not selected_model_key:
        return
    
    model_id = model_options[selected_model_key]
    st.session_state["selected_model_id"] = model_id
    
    # Status filter
    status_filter = st.selectbox("Status", ["pending", "in_progress", "labeled", "skipped"], index=0)
    
    # Load queue
    with st.spinner("Loading labeling queue..."):
        queue_data = api_get(f"/api/v1/active-learning/models/{model_id}/queue?status={status_filter}")
    
    if not queue_data:
        st.info(f"No {status_filter} items in the queue for this model.")
        return
    
    # Display metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total in Queue", len(queue_data))
    
    with col2:
        pending_count = len([item for item in queue_data if item["status"] == "pending"])
        st.metric("Pending", pending_count)
    
    with col3:
        # Count labeled today (mock for MVP)
        labeled_today = len([item for item in queue_data if item["status"] == "labeled"])
        st.metric("Labeled Today", labeled_today)
    
    with col4:
        high_uncertainty = len([item for item in queue_data if (item.get("uncertainty") or 0) > 0.4])
        st.metric("High Uncertainty", high_uncertainty)
    
    # Display queue table
    if queue_data:
        table_data = []
        for item in queue_data:
            smiles = item.get("compound_smiles", "")
            smiles_display = smiles[:30] + "..." if len(smiles) > 30 else smiles
            
            uncertainty = item.get("uncertainty", 0)
            uncertainty_color = get_uncertainty_color(uncertainty)
            
            table_data.append({
                "Compound": item.get("compound_name", item["compound_id"]),
                "SMILES": smiles_display,
                "Prediction": f"{item.get('prediction', 0):.3f}" if item.get('prediction') else "N/A",
                "Uncertainty": f"{uncertainty:.3f}" if uncertainty else "N/A",
                "Priority": f"{item.get('priority_score', 0):.1f}" if item.get('priority_score') else "N/A",
                "Status": item["status"],
                "ID": item["id"]
            })
        
        df = pd.DataFrame(table_data)
        
        # Display table with click handling
        event = st.dataframe(
            df.drop("ID", axis=1),  # Hide ID column
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Handle row selection
        if event.selection and event.selection.rows:
            selected_idx = event.selection.rows[0]
            selected_item_id = table_data[selected_idx]["ID"]
            st.session_state["selected_item_id"] = selected_item_id
            st.session_state["selected_item_data"] = queue_data[selected_idx]
            st.info(f"Selected compound for labeling. Go to 'Labeling Interface' tab.")


def render_labeling_interface_tab():
    """Render Tab 2: Labeling Interface."""
    st.subheader("Labeling Interface")
    
    # Check if item is selected from Tab 1
    selected_item_id = st.session_state.get("selected_item_id")
    selected_item_data = st.session_state.get("selected_item_data")
    
    if not selected_item_id:
        st.info("No compound selected. Please select a compound from the Label Queue tab.")
        return
    
    # Display compound information
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Compound Information")
        st.text(f"Compound: {selected_item_data.get('compound_name', 'N/A')}")
        smiles = selected_item_data.get('compound_smiles', '')
        st.text(f"SMILES: {smiles}")
        
        # Try to display structure (placeholder for now)
        if smiles:
            st.markdown("**Structure:**")
            st.code(smiles, language="text")  # Placeholder - would use mol_to_svg in production
        
        # Properties (would fetch from compound API in production)
        st.markdown("**Properties:**")
        st.text("MW: ~300 Da (estimated)")
        st.text("LogP: ~2.5 (estimated)")
    
    with col2:
        st.markdown("#### Model Prediction")
        prediction = selected_item_data.get('prediction', 0)
        uncertainty = selected_item_data.get('uncertainty', 0)
        
        # Prediction with uncertainty visualization
        st.metric("Prediction", f"{prediction:.3f}")
        st.metric("Uncertainty", f"{uncertainty:.3f}")
        
        # Confidence visualization
        confidence_pct = max(0, min(100, (1 - uncertainty) * 100))
        st.progress(confidence_pct / 100, text=f"Model is {confidence_pct:.0f}% confident")
        
        # Strategy info
        st.text(f"Selection Strategy: {selected_item_data.get('selection_strategy', 'N/A')}")
        st.text(f"Priority Score: {selected_item_data.get('priority_score', 'N/A')}")
    
    # Labeling form
    st.divider()
    st.markdown("#### Submit Label")
    
    with st.form("label_form"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            label_value = st.number_input("Label Value", value=0.0, step=0.1, format="%.3f")
        
        with col2:
            label_source = st.selectbox("Source", ["experimental", "literature", "expert_estimate"])
        
        with col3:
            label_confidence = st.radio("Confidence", ["high", "medium", "low"], horizontal=True)
        
        notes = st.text_area("Notes (optional)", placeholder="Additional context or observations...")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            submit_label = st.form_submit_button("Save Label", type="primary")
        
        with col2:
            skip_item = st.form_submit_button("Skip Item")
        
        with col3:
            next_item = st.form_submit_button("Next Item")
        
        if submit_label:
            payload = {
                "label": label_value,
                "source": label_source,
                "confidence": label_confidence,
                "notes": notes if notes else None
            }
            
            with st.spinner("Submitting label..."):
                result = api_post(f"/api/v1/active-learning/queue/{selected_item_id}/label", payload)
                if result:
                    st.success(f"Label submitted: {label_value}")
                    # Clear selection
                    del st.session_state["selected_item_id"]
                    del st.session_state["selected_item_data"]
                    st.rerun()
        
        if skip_item:
            with st.spinner("Skipping item..."):
                result = api_post(f"/api/v1/active-learning/queue/{selected_item_id}/skip", {})
                if result:
                    st.success("Item skipped")
                    # Clear selection
                    del st.session_state["selected_item_id"]
                    del st.session_state["selected_item_data"]
                    st.rerun()
        
        if next_item:
            # Clear current selection to force selection of next item
            del st.session_state["selected_item_id"]
            del st.session_state["selected_item_data"]
            st.info("Selection cleared. Please select the next compound from the Label Queue tab.")
            st.rerun()


def render_cycle_history_tab():
    """Render Tab 3: Cycle History."""
    st.subheader("Cycle History")
    
    # Model selection
    models = get_ml_models()
    if not models:
        st.info("No ML models available.")
        return
    
    model_options = {f"{m['name']} ({m['model_type']})": m["id"] for m in models}
    selected_model_key = st.selectbox("Select Model for Cycle History", list(model_options.keys()))
    
    if not selected_model_key:
        return
    
    model_id = model_options[selected_model_key]
    
    # Load cycles and stats
    col1, col2 = st.columns(2)
    
    with col1:
        with st.spinner("Loading cycle history..."):
            cycles = api_get(f"/api/v1/active-learning/models/{model_id}/cycles")
    
    with col2:
        with st.spinner("Loading statistics..."):
            stats = api_get(f"/api/v1/active-learning/models/{model_id}/stats")
    
    # Display statistics
    if stats:
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Cycles", stats.get("total_cycles", 0))
        with col2:
            st.metric("Pending Items", stats.get("pending_items", 0))
        with col3:
            st.metric("Labeled Items", stats.get("labeled_items", 0))
        with col4:
            completion_rate = 0
            if stats.get("labeled_items", 0) + stats.get("pending_items", 0) > 0:
                completion_rate = stats.get("labeled_items", 0) / (stats.get("labeled_items", 0) + stats.get("pending_items", 0)) * 100
            st.metric("Completion Rate", f"{completion_rate:.1f}%")
    
    # Display cycles table
    if cycles:
        st.markdown("#### Active Learning Cycles")
        cycle_data = []
        for cycle in cycles:
            cycle_data.append({
                "Cycle #": cycle["cycle_number"],
                "Strategy": cycle.get("selection_strategy", "N/A"),
                "Batch Size": cycle.get("batch_size", 0),
                "Status": cycle["status"],
                "Selected": cycle.get("items_selected", 0),
                "Labeled": cycle.get("items_labeled", 0),
                "Skipped": cycle.get("items_skipped", 0),
                "Started": cycle["started_at"][:10] if cycle.get("started_at") else "N/A"
            })
        
        cycle_df = pd.DataFrame(cycle_data)
        st.dataframe(cycle_df, hide_index=True, use_container_width=True)
        
        # Learning curve visualization
        if len(cycles) > 1:
            st.markdown("#### Learning Curve")
            
            # Extract metrics for plotting (mock data for MVP)
            cycle_numbers = [c["cycle_number"] for c in cycles if c.get("metrics_after")]
            if cycle_numbers:
                # Mock learning curve - in production, extract from metrics_after
                mock_auc = [0.65 + 0.05 * i + 0.02 * (i % 3) for i in range(len(cycle_numbers))]
                
                fig = px.line(
                    x=cycle_numbers,
                    y=mock_auc,
                    title="Model Performance Over Cycles",
                    labels={"x": "Cycle Number", "y": "AUC"}
                )
                fig.update_traces(mode="lines+markers")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No completed cycles with metrics available yet.")
    else:
        st.info("No active learning cycles found for this model.")
    
    # Create new cycle form
    st.divider()
    st.markdown("#### Create New Cycle")
    
    with st.form("create_cycle_form"):
        col1, col2 = st.columns(2)
        
        with col1:
            strategy = st.selectbox(
                "Selection Strategy",
                ["uncertainty", "margin", "entropy", "hybrid", "random"],
                help="Strategy for selecting informative samples"
            )
        
        with col2:
            batch_size = st.slider("Batch Size", min_value=10, max_value=200, value=50, step=10)
        
        submitted = st.form_submit_button("Create Cycle", type="primary")
        
        if submitted:
            payload = {
                "selection_strategy": strategy,
                "batch_size": batch_size
            }
            
            with st.spinner("Creating new cycle..."):
                result = api_post(f"/api/v1/active-learning/models/{model_id}/cycles", payload)
                if result:
                    st.success(f"Created cycle {result.get('cycle_number')} with {strategy} strategy")
                    st.rerun()


def render_sample_selection_tab():
    """Render Tab 4: Sample Selection."""
    st.subheader("Sample Selection")
    
    # Model selection
    models = get_ml_models()
    if not models:
        st.info("No ML models available.")
        return
    
    model_options = {f"{m['name']} ({m['model_type']})": m["id"] for m in models}
    selected_model_key = st.selectbox("Select Model for Sample Selection", list(model_options.keys()))
    
    if not selected_model_key:
        return
    
    model_id = model_options[selected_model_key]
    
    # Pool configuration
    st.markdown("#### Compound Pool Configuration")
    
    pool_type = st.radio(
        "Compound Pool",
        ["All unlabeled compounds", "Specific dataset", "Custom list"],
        help="Choose the pool of compounds to select from"
    )
    
    compound_ids = []
    
    if pool_type == "All unlabeled compounds":
        # Mock compound IDs for MVP
        compound_ids = [str(uuid.uuid4()) for _ in range(1000)]
        st.info(f"Using pool of {len(compound_ids)} unlabeled compounds")
    
    elif pool_type == "Specific dataset":
        # Would query datasets API in production
        dataset_options = ["Dataset A (500 compounds)", "Dataset B (300 compounds)", "Dataset C (200 compounds)"]
        selected_dataset = st.selectbox("Select Dataset", dataset_options)
        if selected_dataset:
            # Mock compound IDs from dataset
            compound_count = int(selected_dataset.split("(")[1].split(" ")[0])
            compound_ids = [str(uuid.uuid4()) for _ in range(compound_count)]
            st.info(f"Using {len(compound_ids)} compounds from {selected_dataset}")
    
    elif pool_type == "Custom list":
        compound_input = st.text_area(
            "Compound IDs (one per line)",
            placeholder="Enter compound UUIDs, one per line"
        )
        if compound_input:
            compound_ids = [line.strip() for line in compound_input.split("\n") if line.strip()]
            st.info(f"Using {len(compound_ids)} custom compounds")
    
    if not compound_ids:
        st.warning("No compounds in pool. Please configure the compound pool above.")
        return
    
    # Selection parameters
    st.markdown("#### Selection Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        strategy = st.selectbox(
            "Selection Strategy",
            ["uncertainty", "margin", "entropy", "hybrid", "random"],
            help="Algorithm for selecting informative samples"
        )
    
    with col2:
        batch_size = st.slider("Batch Size", min_value=10, max_value=min(200, len(compound_ids)), value=50, step=10)
    
    # Strategy explanation
    strategy_explanations = {
        "uncertainty": "Selects compounds with highest prediction uncertainty (ensemble std deviation)",
        "margin": "Selects compounds closest to decision boundary (prediction ≈ 0.5)",
        "entropy": "Selects compounds with highest prediction entropy",
        "hybrid": "Combines uncertainty, margin, and applicability with weighted scoring",
        "random": "Random selection for baseline comparison"
    }
    
    st.info(f"**{strategy.title()} Strategy**: {strategy_explanations[strategy]}")
    
    # Preview
    st.markdown("#### Selection Preview")
    st.text(f"Will select {batch_size} samples from pool of {len(compound_ids)} compounds")
    
    # Selection button
    if st.button("Select Samples", type="primary"):
        payload = {
            "compound_ids": compound_ids,
            "strategy": strategy,
            "batch_size": batch_size
        }
        
        with st.spinner(f"Selecting {batch_size} samples using {strategy} strategy..."):
            result = api_post(f"/api/v1/active-learning/models/{model_id}/select", payload)
            
            if result:
                st.success(f"Selected {len(result)} samples for labeling!")
                
                # Show selection results
                st.markdown("#### Selection Results")
                result_data = []
                for item in result[:10]:  # Show first 10
                    result_data.append({
                        "Compound": item["compound_id"],
                        "Prediction": f"{item.get('prediction', 0):.3f}",
                        "Uncertainty": f"{item.get('uncertainty', 0):.3f}",
                        "Priority": f"{item.get('priority_score', 0):.1f}",
                        "Strategy": item.get("selection_strategy", "N/A")
                    })
                
                if result_data:
                    result_df = pd.DataFrame(result_data)
                    st.dataframe(result_df, hide_index=True, use_container_width=True)
                
                st.info("Samples added to labeling queue. Go to 'Label Queue' tab to start labeling.")
            else:
                st.error("Failed to select samples")


def main():
    """Main Active Learning page."""
    st.title("🔄 Active Learning")
    st.markdown("Iterative model improvement through strategic sample selection and labeling.")
    
    # Create tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "Label Queue", 
        "Labeling Interface", 
        "Cycle History", 
        "Sample Selection"
    ])
    
    with tab1:
        render_label_queue_tab()
    
    with tab2:
        render_labeling_interface_tab()
    
    with tab3:
        render_cycle_history_tab()
    
    with tab4:
        render_sample_selection_tab()


if __name__ == "__main__":
    main()
