"""ML Predictors (Assay Outcome Prediction) page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
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


def render_predictors_page() -> None:
    """Render the ML Predictors page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🤖 ML Predictors")
    st.caption("Train and run assay outcome prediction models for screening prioritization.")
    
    tab1, tab2, tab3 = st.tabs(["Models", "Train Model", "Run Predictions"])
    
    with tab1:
        render_models_tab()
    
    with tab2:
        render_train_tab()
    
    with tab3:
        render_predict_tab()


def render_models_tab() -> None:
    """Render the models list tab."""
    st.subheader("Available Predictors")
    
    # Optional program filter
    program_id = st.text_input(
        "Filter by Program ID (optional)",
        placeholder="UUID of program...",
        help="Leave empty to show all models",
    )
    
    if st.button("Load Models"):
        try:
            params = f"?program_id={program_id}" if program_id else ""
            response = _api_get(f"/api/v1/predictors{params}")
            st.session_state["predictors_models"] = response
        except Exception as e:
            st.error(f"Failed to load models: {e}")
            return
    
    models_data = st.session_state.get("predictors_models")
    
    if not models_data:
        st.info("Click 'Load Models' to view available predictors.")
        return
    
    models = models_data.get("models", [])
    
    if not models:
        st.info("No predictors found.")
        return
    
    # Display models table
    model_rows = []
    for m in models:
        perf = m.get("performance", {})
        model_rows.append({
            "Model ID": str(m.get("model_id", ""))[:16] + "...",
            "Name": m.get("name", ""),
            "Assay Type": m.get("assay_type", ""),
            "Version": m.get("version", ""),
            "Accuracy": f"{perf.get('accuracy', 0.0):.3f}" if perf.get('accuracy') else "N/A",
            "ROC-AUC": f"{perf.get('roc_auc', 0.0):.3f}" if perf.get('roc_auc') else "N/A",
            "Created": m.get("created_at", "")[:10] if m.get("created_at") else "",
        })
    
    df = pd.DataFrame(model_rows)
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    st.caption(f"Total models: {models_data.get('total_models', 0)}")


def render_train_tab() -> None:
    """Render the model training tab."""
    st.subheader("Train New Predictor")
    
    st.markdown("**Note**: Training requires HTS campaign data with compound activities.")
    
    program_id = st.text_input(
        "Program ID",
        placeholder="UUID of program with screening data...",
        help="Program must have HTS campaigns or assay results",
    )
    
    assay_type = st.text_input(
        "Assay Type",
        placeholder="e.g., cell_viability, binding_assay...",
        help="Must match assay type in HTS campaigns",
    )
    
    col1, col2 = st.columns(2)
    with col1:
        min_actives = st.number_input(
            "Min Active Compounds",
            min_value=10,
            max_value=10000,
            value=50,
            step=10,
            help="Minimum number of active compounds needed",
        )
    
    with col2:
        min_inactives = st.number_input(
            "Min Inactive Compounds",
            min_value=10,
            max_value=10000,
            value=50,
            step=10,
            help="Minimum number of inactive compounds needed",
        )
    
    features = st.multiselect(
        "Features",
        options=["Morgan", "MACCS", "RDKit2D"],
        default=["Morgan"],
        help="Molecular fingerprint types to use",
    )
    
    if st.button("Train Model", type="primary", disabled=not (program_id and assay_type)):
        with st.spinner("Training model... This may take several minutes."):
            try:
                result = _api_post(
                    "/api/v1/predictors/train",
                    {
                        "program_id": program_id,
                        "assay_type": assay_type,
                        "features": features or ["Morgan"],
                        "min_actives": min_actives,
                        "min_inactives": min_inactives,
                    },
                    timeout=300,
                )
                st.session_state["training_result"] = result
            except Exception as e:
                st.error(f"Training failed: {e}")
                return
    
    # Display training result
    result = st.session_state.get("training_result")
    if result:
        st.success(f"Model trained! ID: {result.get('model_id')}")
        
        perf = result.get("model_performance", {})
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Accuracy", f"{perf.get('accuracy', 0.0):.3f}")
        with col2:
            st.metric("ROC-AUC", f"{perf.get('roc_auc', 0.0):.3f}")
        with col3:
            st.metric("Training Time", f"{result.get('training_time_seconds', 0.0):.1f}s")
        
        stats = result.get("training_stats", {})
        st.markdown("**Training Data**:")
        st.markdown(f"- Total compounds: {stats.get('total_compounds', 0)}")
        st.markdown(f"- Active: {stats.get('active_compounds', 0)}")
        st.markdown(f"- Inactive: {stats.get('inactive_compounds', 0)}")


def render_predict_tab() -> None:
    """Render the predictions tab."""
    st.subheader("Run Predictions")
    
    model_id = st.text_input(
        "Model ID",
        placeholder="UUID of trained model...",
        help="Get model ID from the Models tab",
    )
    
    smiles_input = st.text_area(
        "Compound SMILES (one per line)",
        height=200,
        placeholder="CCO\nCC(C)O\nC1CCCCC1",
        help="Enter SMILES strings for compounds to predict",
    )
    
    if st.button("Predict", type="primary", disabled=not (model_id and smiles_input)):
        smiles_list = [s.strip() for s in smiles_input.strip().split("\n") if s.strip()]
        
        with st.spinner(f"Running predictions for {len(smiles_list)} compounds..."):
            try:
                result = _api_post(
                    f"/api/v1/predictors/{model_id}/predict",
                    {"compound_smiles": smiles_list},
                    timeout=120,
                )
                st.session_state["prediction_result"] = result
            except Exception as e:
                st.error(f"Prediction failed: {e}")
                return
    
    # Display prediction result
    result = st.session_state.get("prediction_result")
    if result:
        predictions = result.get("predictions", [])
        
        if not predictions:
            st.warning("No predictions returned.")
            return
        
        st.success(f"Predicted {result.get('total_predictions', 0)} compounds in {result.get('processing_time_seconds', 0.0):.2f}s")
        
        # Display predictions table
        pred_rows = []
        for p in predictions:
            pred_rows.append({
                "SMILES": p.get("compound_smiles", ""),
                "Prediction": "Active" if p.get("prediction") else "Inactive",
                "Probability Active": f"{p.get('probability_active', 0.0):.3f}",
                "Confidence": f"{p.get('confidence', 0.0):.3f}",
            })
        
        df = pd.DataFrame(pred_rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        # Download results
        csv = df.to_csv(index=False)
        st.download_button(
            "📥 Download Predictions (CSV)",
            data=csv,
            file_name="predictions.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    render_predictors_page()

