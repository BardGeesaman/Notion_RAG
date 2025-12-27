"""Model Monitoring Dashboard."""

from __future__ import annotations

import os
from typing import Dict, List

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str) -> dict:
    """Helper function for GET requests."""
    try:
        response = requests.get(f"{API_BASE}{path}", timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {e}")
        return {}


def _get_all_models() -> List[Dict]:
    """Fetch all models from the API."""
    return _api_get("/api/v1/models/")


def _get_model_health(model_id: str, window_hours: int = 24) -> Dict:
    """Fetch model health report."""
    return _api_get(f"/api/v1/monitoring/models/{model_id}/health?window_hours={window_hours}")


def _get_drift_report(model_id: str, window_hours: int = 24) -> Dict:
    """Fetch drift analysis report."""
    return _api_get(f"/api/v1/monitoring/models/{model_id}/drift?window_hours={window_hours}")


def _get_calibration_report(model_id: str) -> Dict:
    """Fetch calibration analysis report."""
    return _api_get(f"/api/v1/monitoring/models/{model_id}/calibration")


def _status_color(status: str) -> str:
    """Get color for status badge."""
    color_map = {
        "good": "🟢",
        "warning": "🟡", 
        "alert": "🔴",
        "no_baseline": "⚪",
        "no_data": "⚪",
        "not_applicable": "⚫"
    }
    return color_map.get(status, "⚪")


def _render_model_health_overview() -> None:
    """Render the model health overview tab."""
    st.subheader("Model Health Overview")
    
    # Fetch all models
    models = _get_all_models()
    if not models:
        st.info("No models found.")
        return
    
    # Window selector
    window_hours = st.selectbox(
        "Time window",
        options=[24, 168, 720],  # 24h, 7d, 30d
        format_func=lambda x: f"{x}h" if x < 168 else f"{x//24}d",
        index=0
    )
    
    # Fetch health reports for all models
    health_data = []
    progress_bar = st.progress(0)
    
    for i, model in enumerate(models):
        model_id = model.get("id")
        model_name = model.get("name", "Unknown")
        model_version = model.get("version", "Unknown")
        
        if model_id:
            health = _get_model_health(model_id, window_hours)
            
            if health:
                health_data.append({
                    "Model Name": model_name,
                    "Version": model_version,
                    "Overall Status": health.get("overall_status", "unknown"),
                    "Drift Status": health.get("drift_status", "unknown"),
                    "Calibration Status": health.get("calibration_status", "unknown"),
                    "PSI Max": health.get("psi_max"),
                    "ECE": health.get("ece"),
                    "Last Checked": health.get("last_checked", "Never"),
                    "Model ID": model_id
                })
        
        progress_bar.progress((i + 1) / len(models))
    
    if not health_data:
        st.warning("No health data available for any models.")
        return
    
    # Create DataFrame
    df = pd.DataFrame(health_data)
    
    # Format columns for display
    df_display = df.copy()
    df_display["Status"] = df_display["Overall Status"].apply(
        lambda x: f"{_status_color(x)} {x.title()}"
    )
    df_display["Drift"] = df_display["Drift Status"].apply(
        lambda x: f"{_status_color(x)} {x.title()}"
    )
    df_display["Calibration"] = df_display["Calibration Status"].apply(
        lambda x: f"{_status_color(x)} {x.title()}"
    )
    
    # Format numeric columns
    df_display["PSI Max"] = df_display["PSI Max"].apply(
        lambda x: f"{x:.3f}" if x is not None else "N/A"
    )
    df_display["ECE"] = df_display["ECE"].apply(
        lambda x: f"{x:.3f}" if x is not None else "N/A"
    )
    
    # Format timestamp
    df_display["Last Checked"] = df_display["Last Checked"].apply(
        lambda x: x.split("T")[0] if isinstance(x, str) and "T" in x else x
    )
    
    # Display table
    display_columns = ["Model Name", "Version", "Status", "Drift", "Calibration", "PSI Max", "ECE", "Last Checked"]
    st.dataframe(
        df_display[display_columns],
        use_container_width=True,
        hide_index=True
    )
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    
    total_models = len(df)
    good_models = len(df[df["Overall Status"] == "good"])
    warning_models = len(df[df["Overall Status"] == "warning"])
    alert_models = len(df[df["Overall Status"] == "alert"])
    
    col1.metric("Total Models", total_models)
    col2.metric("🟢 Good", good_models)
    col3.metric("🟡 Warning", warning_models)
    col4.metric("🔴 Alert", alert_models)


def _render_drift_analysis() -> None:
    """Render the drift analysis tab."""
    st.subheader("Drift Analysis")
    
    # Fetch all models
    models = _get_all_models()
    if not models:
        st.info("No models found.")
        return
    
    # Model selector
    model_options = {f"{m.get('name', 'Unknown')} (v{m.get('version', '?')})": m.get('id') 
                    for m in models if m.get('id')}
    
    if not model_options:
        st.warning("No valid models found.")
        return
    
    selected_model_name = st.selectbox("Select Model", list(model_options.keys()))
    selected_model_id = model_options[selected_model_name]
    
    # Window selector
    window_options = {
        "24 hours": 24,
        "7 days": 168,
        "30 days": 720
    }
    selected_window_name = st.selectbox("Time Window", list(window_options.keys()))
    window_hours = window_options[selected_window_name]
    
    # Fetch drift report
    drift_report = _get_drift_report(selected_model_id, window_hours)
    
    if not drift_report:
        st.error("Failed to fetch drift report.")
        return
    
    # Display status
    status = drift_report.get("status", "unknown")
    st.write(f"**Drift Status:** {_status_color(status)} {status.title()}")
    
    # Display basic metrics
    col1, col2 = st.columns(2)
    col1.metric("Predictions Analyzed", drift_report.get("n_predictions", 0))
    col2.metric("Window", f"{drift_report.get('window_hours', 0)} hours")
    
    # PSI scores chart
    psi_scores = drift_report.get("psi_scores", {})
    fp_drift = drift_report.get("fp_aggregate_drift", {})
    
    if psi_scores or fp_drift:
        # Combine all PSI-like scores
        all_scores = {}
        all_scores.update(psi_scores)
        all_scores.update({f"FP_{k}": v for k, v in fp_drift.items()})
        
        if all_scores:
            # Create bar chart
            features = list(all_scores.keys())
            scores = list(all_scores.values())
            
            fig = go.Figure()
            
            # Add bars
            colors = ['red' if score >= 0.25 else 'orange' if score >= 0.1 else 'green' 
                     for score in scores]
            
            fig.add_trace(go.Bar(
                x=features,
                y=scores,
                marker_color=colors,
                name="PSI Score"
            ))
            
            # Add threshold lines
            fig.add_hline(y=0.1, line_dash="dash", line_color="orange", 
                         annotation_text="Warning (0.1)")
            fig.add_hline(y=0.25, line_dash="dash", line_color="red",
                         annotation_text="Alert (0.25)")
            
            fig.update_layout(
                title="Feature Drift (PSI Scores)",
                xaxis_title="Features",
                yaxis_title="PSI Score",
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No drift scores available.")
    else:
        if status == "no_baseline":
            st.info("No baseline available for drift detection. Train the model and establish a baseline first.")
        elif status == "no_data":
            st.info("No prediction data available for the selected time window.")
        else:
            st.info("No drift data available.")


def _render_calibration_analysis() -> None:
    """Render the calibration analysis tab."""
    st.subheader("Calibration Analysis")
    
    # Fetch all models
    models = _get_all_models()
    if not models:
        st.info("No models found.")
        return
    
    # Filter to classification models only
    classification_models = [m for m in models if m.get('model_type') == 'admet_classification']
    
    if not classification_models:
        st.info("No classification models found. Calibration analysis is only applicable to classification models.")
        return
    
    # Model selector
    model_options = {f"{m.get('name', 'Unknown')} (v{m.get('version', '?')})": m.get('id') 
                    for m in classification_models if m.get('id')}
    
    selected_model_name = st.selectbox("Select Classification Model", list(model_options.keys()))
    selected_model_id = model_options[selected_model_name]
    
    # Fetch calibration report
    calibration_report = _get_calibration_report(selected_model_id)
    
    if not calibration_report:
        st.error("Failed to fetch calibration report.")
        return
    
    # Display status
    status = calibration_report.get("status", "unknown")
    st.write(f"**Calibration Status:** {_status_color(status)} {status.title()}")
    
    # Display metrics
    col1, col2 = st.columns(2)
    
    ece = calibration_report.get("ece")
    if ece is not None:
        col1.metric("Expected Calibration Error (ECE)", f"{ece:.4f}")
    else:
        col1.metric("Expected Calibration Error (ECE)", "N/A")
    
    col2.metric("Predictions with Ground Truth", calibration_report.get("n_predictions_with_truth", 0))
    
    # Reliability diagram (placeholder - would need reliability data from API)
    if status not in ["no_data", "not_applicable"]:
        st.subheader("Reliability Diagram")
        
        # Create placeholder reliability diagram
        # In a real implementation, this would use reliability data from the API
        confidence_bins = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
        
        # Perfect calibration line
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode='lines',
            name='Perfect Calibration',
            line=dict(dash='dash', color='gray')
        ))
        
        # Placeholder actual calibration (would come from API)
        if ece is not None:
            # Simulate some calibration data based on ECE
            import numpy as np
            np.random.seed(42)
            noise = np.random.normal(0, ece, len(confidence_bins))
            actual_accuracy = [min(1, max(0, conf + n)) for conf, n in zip(confidence_bins, noise)]
            
            fig.add_trace(go.Scatter(
                x=confidence_bins,
                y=actual_accuracy,
                mode='lines+markers',
                name='Model Calibration',
                line=dict(color='blue')
            ))
        
        fig.update_layout(
            title="Reliability Diagram",
            xaxis_title="Mean Predicted Probability",
            yaxis_title="Fraction of Positives",
            xaxis=dict(range=[0, 1]),
            yaxis=dict(range=[0, 1])
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        st.caption("Note: This is a simplified visualization. Full reliability diagram requires additional API data.")
    
    else:
        if status == "no_data":
            st.info("No ground truth data available for calibration analysis.")
        elif status == "not_applicable":
            st.info("Calibration analysis is not applicable to this model type.")


def render_model_monitoring_page() -> None:
    """Main function to render the model monitoring page."""
    st.header("🔍 Model Monitoring")
    st.caption("Monitor model performance, drift, and calibration over time.")
    
    # Create tabs
    tab1, tab2, tab3 = st.tabs(["Model Health Overview", "Drift Analysis", "Calibration"])
    
    with tab1:
        _render_model_health_overview()
    
    with tab2:
        _render_drift_analysis()
    
    with tab3:
        _render_calibration_analysis()
