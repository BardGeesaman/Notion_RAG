"""Deep Learning Toxicity Predictions Dashboard."""
import os
from typing import List, Dict, Any

import streamlit as st
import requests
import pandas as pd
import plotly.express as px

API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def get_auth_headers() -> Dict[str, str]:
    """Get authentication headers."""
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


def api_post(path: str, payload: Dict[str, Any], params: Dict[str, Any] = None) -> Dict[str, Any]:
    """Make authenticated POST request to API."""
    try:
        response = requests.post(
            f"{API_BASE}{path}", 
            json=payload, 
            params=params,
            headers=get_auth_headers(), 
            timeout=60.0
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        st.error(f"API request failed: {str(e)}")
        return {}


def get_gnn_models() -> List[Dict[str, Any]]:
    """Fetch available GNN models."""
    models_data = api_get("/api/admet/gnn-models")
    return models_data.get("models", [])


def predict_gnn(smiles_list: List[str], endpoints: List[str], with_uncertainty: bool = True) -> Dict[str, Any]:
    """Get GNN predictions."""
    params = {
        "endpoints": endpoints,
        "with_uncertainty": with_uncertainty
    }
    return api_post("/api/admet/predict-gnn", {"smiles_list": smiles_list}, params)


def compare_models(smiles_list: List[str], endpoint: str) -> Dict[str, Any]:
    """Compare XGBoost vs GNN."""
    params = {"endpoint": endpoint}
    return api_post("/api/admet/compare", {"smiles_list": smiles_list}, params)


def get_risk_color(value: float) -> str:
    """Get risk color for toxicity score."""
    if value > 0.7:
        return "🔴"  # High risk
    elif value > 0.4:
        return "🟡"  # Moderate risk
    else:
        return "🟢"  # Low risk


def render_predict_tab():
    """Render Tab 1: Predict."""
    st.subheader("🔮 GNN Toxicity Prediction")
    st.markdown("Graph Neural Network models for molecular toxicity screening with uncertainty quantification.")
    
    # Input section
    col1, col2 = st.columns([2, 1])
    
    with col1:
        smiles_input = st.text_area(
            "Enter SMILES (one per line)",
            value="CCO\nc1ccccc1\nCC(=O)Oc1ccccc1C(=O)O\nCCN(CC)CC",
            height=150,
            help="Enter molecular SMILES strings, one per line"
        )
    
    with col2:
        # Get available models
        models = get_gnn_models()
        available = [m["endpoint"] for m in models if m.get("available")]
        
        if available:
            selected_endpoints = st.multiselect(
                "Select Toxicity Endpoints",
                options=available,
                default=available[:1],
                help="Choose which toxicity endpoints to predict"
            )
        else:
            st.warning("⚠️ No GNN models available. Train models first.")
            selected_endpoints = []
        
        with_uncertainty = st.checkbox(
            "Include uncertainty estimation", 
            value=True,
            help="Use MC Dropout for uncertainty quantification"
        )
        
        n_samples = st.slider(
            "MC Dropout samples", 
            min_value=5, 
            max_value=20, 
            value=10,
            help="More samples = better uncertainty, slower inference"
        ) if with_uncertainty else 10
    
    if st.button("🔮 Predict Toxicity", type="primary", disabled=not selected_endpoints):
        smiles_list = [s.strip() for s in smiles_input.strip().split("\n") if s.strip()]
        
        if not smiles_list:
            st.error("Please enter at least one SMILES string")
            return
        
        with st.spinner(f"Running GNN inference on {len(smiles_list)} molecules..."):
            results = predict_gnn(smiles_list, selected_endpoints, with_uncertainty)
        
        if results and results.get("predictions"):
            st.success(f"✅ Predicted {len(results['predictions'])} molecules")
            
            # Summary metrics
            col1, col2, col3 = st.columns(3)
            predictions = results["predictions"]
            
            with col1:
                st.metric("Molecules Processed", len(predictions))
            
            with col2:
                error_count = sum(1 for p in predictions 
                                for ep in p.get("endpoints", {}).values() 
                                if "error" in ep)
                st.metric("Prediction Errors", error_count)
            
            with col3:
                if selected_endpoints and with_uncertainty:
                    avg_uncertainty = 0
                    count = 0
                    for p in predictions:
                        for ep in selected_endpoints:
                            ep_result = p.get("endpoints", {}).get(ep, {})
                            if "uncertainty" in ep_result:
                                avg_uncertainty += ep_result["uncertainty"]
                                count += 1
                    if count > 0:
                        avg_uncertainty /= count
                        st.metric("Avg Uncertainty", f"±{avg_uncertainty:.3f}")
            
            # Display results
            for pred in predictions:
                smiles = pred["smiles"]
                with st.expander(f"📝 {smiles[:50]}{'...' if len(smiles) > 50 else ''}", expanded=True):
                    cols = st.columns(len(selected_endpoints))
                    
                    for i, endpoint in enumerate(selected_endpoints):
                        with cols[i]:
                            ep_result = pred.get("endpoints", {}).get(endpoint, {})
                            
                            if "error" in ep_result:
                                st.error(f"❌ {ep_result['error']}")
                            else:
                                value = ep_result.get("prediction", 0)
                                unc = ep_result.get("uncertainty")
                                label = ep_result.get("label", "")
                                
                                # Color-coded metric with risk assessment
                                risk_icon = get_risk_color(value)
                                
                                if value > 0.7:
                                    delta = "⚠️ High Risk"
                                    delta_color = "inverse"
                                elif value > 0.4:
                                    delta = "⚡ Moderate"
                                    delta_color = "normal"
                                else:
                                    delta = "✅ Low Risk"
                                    delta_color = "normal"
                                
                                st.metric(
                                    f"{risk_icon} {endpoint.upper()}", 
                                    f"{value:.1%}", 
                                    delta=delta,
                                    delta_color=delta_color
                                )
                                
                                if unc is not None:
                                    st.caption(f"Uncertainty: ±{unc:.3f}")
                                if label:
                                    st.caption(f"Classification: {label}")


def render_compare_tab():
    """Render Tab 2: Compare Models."""
    st.subheader("📊 XGBoost vs GNN Comparison")
    st.markdown("Compare traditional machine learning (XGBoost ensemble) with deep learning (Graph Neural Network) approaches.")
    
    col1, col2 = st.columns(2)
    
    with col1:
        compare_smiles = st.text_input(
            "SMILES for comparison",
            value="c1ccc(cc1)N",
            help="Enter a single SMILES string to compare both models"
        )
    
    with col2:
        models = get_gnn_models()
        available = [m["endpoint"] for m in models if m.get("available")]
        
        if available:
            compare_endpoint = st.selectbox(
                "Toxicity Endpoint", 
                options=available,
                help="Select endpoint for model comparison"
            )
        else:
            st.warning("No GNN models available for comparison")
            return
    
    if st.button("📊 Compare Models", type="primary"):
        if not compare_smiles.strip():
            st.error("Please enter a SMILES string")
            return
        
        with st.spinner("Running model comparison..."):
            results = compare_models([compare_smiles.strip()], compare_endpoint)
        
        if results and results.get("comparisons"):
            comp = results["comparisons"][0]
            
            st.markdown("### Model Comparison Results")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("#### 🌲 XGBoost Ensemble")
                xgb = comp.get("xgboost", {})
                if "error" in xgb:
                    st.warning(f"❌ {xgb['error']}")
                else:
                    prob = xgb.get("probability", xgb.get("value", 0))
                    risk_icon = get_risk_color(prob)
                    st.metric("Prediction", f"{risk_icon} {prob:.1%}")
                    
                    if "uncertainty" in xgb:
                        st.caption(f"Uncertainty: ±{xgb['uncertainty']:.3f}")
                    if "applicability" in xgb:
                        st.caption(f"Applicability: {xgb['applicability']:.3f}")
            
            with col2:
                st.markdown("#### 🧬 Graph Neural Network")
                gnn = comp.get("gnn", {})
                if "error" in gnn:
                    st.warning(f"❌ {gnn['error']}")
                else:
                    pred = gnn.get("prediction", 0)
                    risk_icon = get_risk_color(pred)
                    st.metric("Prediction", f"{risk_icon} {pred:.1%}")
                    
                    if "uncertainty" in gnn:
                        st.caption(f"Uncertainty: ±{gnn['uncertainty']:.3f}")
                    if "confidence" in gnn:
                        st.caption(f"Confidence: {gnn['confidence']:.3f}")
                    if "label" in gnn:
                        st.caption(f"Classification: {gnn['label']}")
            
            # Model comparison insights
            if "error" not in xgb and "error" not in gnn:
                st.markdown("### 🔍 Model Insights")
                
                xgb_pred = xgb.get("probability", xgb.get("value", 0))
                gnn_pred = gnn.get("prediction", 0)
                diff = abs(xgb_pred - gnn_pred)
                
                if diff < 0.1:
                    st.success(f"✅ **Models agree** (difference: {diff:.1%})")
                elif diff < 0.3:
                    st.warning(f"⚡ **Moderate disagreement** (difference: {diff:.1%})")
                else:
                    st.error(f"⚠️ **Strong disagreement** (difference: {diff:.1%}) - requires expert review")


def render_model_info_tab():
    """Render Tab 3: Model Info."""
    st.subheader("ℹ️ GNN Model Information")
    st.markdown("Status and performance metrics for available Graph Neural Network models.")
    
    models = get_gnn_models()
    
    if not models:
        st.info("🔍 No models found. Train models using the command below:")
        st.code("python scripts/train_gnn_models.py --endpoint herg --epochs 50")
        return
    
    # Summary metrics
    available_count = sum(1 for m in models if m.get("available"))
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Total Endpoints", len(models))
    with col2:
        st.metric("Available Models", available_count)
    with col3:
        st.metric("Training Required", len(models) - available_count)
    
    # Model details
    for model in models:
        endpoint = model["endpoint"]
        available = model.get("available", False)
        
        with st.expander(f"{'✅' if available else '❌'} {endpoint.upper()} - {_get_endpoint_description(endpoint)}", expanded=available):
            if available:
                metrics = model.get("metrics", {})
                
                # Performance metrics
                if metrics:
                    st.markdown("**Performance Metrics:**")
                    cols = st.columns(4)
                    
                    if "auc" in metrics:
                        cols[0].metric("AUC", f"{metrics['auc']:.1%}")
                    if "accuracy" in metrics:
                        cols[1].metric("Accuracy", f"{metrics['accuracy']:.1%}")
                    if "rmse" in metrics:
                        cols[2].metric("RMSE", f"{metrics['rmse']:.3f}")
                    if "r2" in metrics:
                        cols[3].metric("R²", f"{metrics['r2']:.3f}")
                
                # Training info
                if model.get("trained_at"):
                    st.caption(f"🕒 Trained: {model['trained_at'][:19].replace('T', ' ')}")
                
                task_type = model.get("task_type", "unknown")
                st.caption(f"📊 Task Type: {task_type}")
                
                # Model details button
                if st.button(f"📋 View Details", key=f"details_{endpoint}"):
                    model_info = api_get(f"/api/admet/gnn-models/{endpoint}/info")
                    if model_info:
                        st.json(model_info)
            else:
                error_msg = model.get("error", "Model not trained")
                st.warning(f"⚠️ {error_msg}")
                
                # Training command
                st.markdown("**Training Command:**")
                st.code(f"python scripts/train_gnn_models.py --endpoint {endpoint} --epochs 50")


def _get_endpoint_description(endpoint: str) -> str:
    """Get description for toxicity endpoint."""
    descriptions = {
        "herg": "hERG channel inhibition (cardiotoxicity)",
        "ames": "Mutagenicity (Ames test)",
        "dili": "Drug-induced liver injury",
        "ld50": "Acute oral toxicity (LD50)",
        "clintox": "Clinical trial toxicity"
    }
    return descriptions.get(endpoint, "Unknown endpoint")


def main():
    """Main GNN Toxicity dashboard."""
    st.title("🧬 Deep Learning Toxicity Predictions")
    st.markdown("Graph Neural Network (GNN) models for molecular toxicity screening with uncertainty quantification.")
    
    # Create tabs
    tab1, tab2, tab3 = st.tabs(["🔮 Predict", "📊 Compare Models", "ℹ️ Model Info"])
    
    with tab1:
        render_predict_tab()
    
    with tab2:
        render_compare_tab()
    
    with tab3:
        render_model_info_tab()


if __name__ == "__main__":
    main()
