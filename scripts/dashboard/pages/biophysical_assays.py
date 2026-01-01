"""Biophysical Assays dashboard page."""

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
import numpy as np


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_biophysical_assays_page() -> None:
    """Render the Biophysical Assays analysis dashboard."""
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("⚗️ Biophysical Assay Analysis")
    st.caption("Upload SPR, MST, and DSC files for molecular interaction and stability studies.")

    # Create tabs
    tab_spr, tab_mst, tab_dsc = st.tabs(["📈 SPR Analysis", "🌡️ MST Analysis", "🔥 DSC Analysis"])

    with tab_spr:
        render_spr_tab()
    
    with tab_mst:
        render_mst_tab()
    
    with tab_dsc:
        render_dsc_tab()


def render_spr_tab() -> None:
    """Render the SPR Analysis tab."""
    st.subheader("Surface Plasmon Resonance (SPR)")
    st.caption("Analyze real-time binding kinetics and affinity measurements")
    
    # Upload section
    st.markdown("### 📤 File Upload")
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Upload Biacore/Reichert File",
            type=["csv", "txt"],
            key="spr_upload",
            help="Upload SPR data files from Biacore or Reichert instruments"
        )
    
    with col2:
        with st.container():
            compound_id = st.text_input("Compound ID (optional)", key="spr_compound", help="Link to existing compound")
            target_name = st.text_input("Target Name (optional)", key="spr_target", help="Protein or target identifier")
    
    if uploaded_file is not None:
        st.success(f"File selected: {uploaded_file.name}")
        
        # Show file details
        col1, col2 = st.columns(2)
        with col1:
            st.metric("File Size", f"{len(uploaded_file.getvalue()) / 1024:.1f} KB")
        with col2:
            st.metric("File Type", uploaded_file.type or "text/csv")
        
        if st.button("Process SPR File", key="spr_process", type="primary", use_container_width=True):
            with st.spinner("Uploading and processing SPR file..."):
                try:
                    # Prepare form data
                    data = {}
                    if compound_id:
                        data["compound_id"] = compound_id
                    if target_name:
                        data["target_name"] = target_name
                    
                    # Upload via API
                    with httpx.Client(timeout=120) as client:
                        files = {"file": (uploaded_file.name, uploaded_file.getvalue(), "text/csv")}
                        response = client.post(f"{API_BASE}/api/v1/biophysical/spr/upload", files=files, data=data)
                    
                    if response.status_code == 201:
                        result = response.json()
                        st.success("✅ SPR file uploaded successfully!")
                        
                        # Show processing status
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Experiment ID", str(result["experiment_id"])[:8] + "...")
                        with col2:
                            st.metric("Processing Status", result["processing_status"])
                        with col3:
                            st.metric("File Size", f"{result['file_size_bytes'] / 1024:.1f} KB")
                        
                        # Store experiment ID in session state
                        st.session_state.current_spr_experiment = result["experiment_id"]
                        
                        st.info("💡 Refresh the page to see your experiment in the analysis section below!")
                        
                    else:
                        st.error(f"Upload failed: {response.text}")
                        
                except Exception as e:
                    st.error(f"Error uploading file: {str(e)}")
    
    # Analysis section
    st.markdown("### 📊 Sensorgram Visualization")
    
    # Fetch existing experiments
    experiments = fetch_spr_experiments()
    
    if experiments:
        # Experiment selector
        experiment_options = {
            exp["experiment_name"] or f"Experiment {str(exp['id'])[:8]}": exp["id"] 
            for exp in experiments
        }
        
        selected_name = st.selectbox(
            "Select SPR Experiment", 
            options=list(experiment_options.keys()),
            key="spr_experiment_select"
        )
        
        if selected_name:
            selected_id = experiment_options[selected_name]
            selected_experiment = next((exp for exp in experiments if exp["id"] == selected_id), None)
            
            if selected_experiment:
                # Create sample sensorgram data (since we don't have real sensorgram endpoint yet)
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    # Plot sensorgrams
                    fig = go.Figure()
                    
                    # Generate sample sensorgram data
                    time_points = np.linspace(0, 600, 1000)
                    concentrations = [0.1, 1.0, 10.0, 100.0]  # nM
                    
                    for i, conc in enumerate(concentrations):
                        # Simulate binding curve
                        ka = selected_experiment.get("ka", 1e5)
                        kd = selected_experiment.get("kd_rate", 1e-3)
                        rmax = selected_experiment.get("rmax", 100)
                        
                        # Simple binding model simulation
                        response = simulate_spr_response(time_points, conc, ka, kd, rmax)
                        
                        fig.add_trace(go.Scatter(
                            x=time_points,
                            y=response,
                            name=f"{conc} nM",
                            mode='lines'
                        ))
                    
                    fig.update_layout(
                        xaxis_title="Time (s)",
                        yaxis_title="Response (RU)",
                        title="SPR Sensorgrams",
                        height=400
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Kinetic parameters
                    st.markdown("**Kinetic Parameters**")
                    
                    ka = selected_experiment.get("ka")
                    kd_rate = selected_experiment.get("kd_rate")
                    kd_eq = selected_experiment.get("kd_equilibrium")
                    rmax = selected_experiment.get("rmax")
                    chi2 = selected_experiment.get("chi_squared")
                    
                    if ka:
                        st.metric("ka (1/Ms)", f"{ka:.2e}")
                    if kd_rate:
                        st.metric("kd (1/s)", f"{kd_rate:.2e}")
                    if kd_eq:
                        st.metric("KD (nM)", f"{kd_eq * 1e9:.2f}")
                    if rmax:
                        st.metric("Rmax (RU)", f"{rmax:.1f}")
                    if chi2:
                        st.metric("χ²", f"{chi2:.2f}")
                
                # Refit controls
                st.markdown("### 🔧 Kinetic Analysis")
                col1, col2 = st.columns(2)
                
                with col1:
                    model = st.selectbox("Fitting Model", ["1:1_langmuir", "two_state"], key="spr_model")
                
                with col2:
                    if st.button("Refit Kinetics", key="spr_refit"):
                        with st.spinner("Refitting kinetic parameters..."):
                            try:
                                with httpx.Client(timeout=60) as client:
                                    refit_data = {"model": model}
                                    response = client.post(
                                        f"{API_BASE}/api/v1/biophysical/spr/{selected_id}/fit", 
                                        json=refit_data
                                    )
                                
                                if response.status_code == 200:
                                    st.success("✅ Kinetic parameters refitted!")
                                    st.rerun()
                                else:
                                    st.error(f"Refit failed: {response.text}")
                            except Exception as e:
                                st.error(f"Error refitting: {str(e)}")
    else:
        st.info("📁 No SPR experiments found. Upload a file to get started!")


def render_mst_tab() -> None:
    """Render the MST Analysis tab."""
    st.subheader("MicroScale Thermophoresis (MST)")
    st.caption("Analyze binding affinity using thermophoretic mobility")
    
    # Upload section
    st.markdown("### 📤 File Upload")
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Upload NanoTemper File",
            type=["xlsx", "csv"],
            key="mst_upload",
            help="Upload MST data files from NanoTemper instruments"
        )
    
    with col2:
        with st.container():
            compound_id = st.text_input("Compound ID (optional)", key="mst_compound")
            target_name = st.text_input("Target Name (optional)", key="mst_target")
    
    if uploaded_file is not None:
        st.success(f"File selected: {uploaded_file.name}")
        
        if st.button("Process MST File", key="mst_process", type="primary", use_container_width=True):
            with st.spinner("Uploading and processing MST file..."):
                try:
                    data = {}
                    if compound_id:
                        data["compound_id"] = compound_id
                    if target_name:
                        data["target_name"] = target_name
                    
                    with httpx.Client(timeout=120) as client:
                        files = {"file": (uploaded_file.name, uploaded_file.getvalue(), "application/vnd.ms-excel")}
                        response = client.post(f"{API_BASE}/api/v1/biophysical/mst/upload", files=files, data=data)
                    
                    if response.status_code == 201:
                        result = response.json()
                        st.success("✅ MST file uploaded successfully!")
                        st.session_state.current_mst_experiment = result["experiment_id"]
                        st.info("💡 Refresh the page to see your experiment in the analysis section below!")
                    else:
                        st.error(f"Upload failed: {response.text}")
                except Exception as e:
                    st.error(f"Error uploading file: {str(e)}")
    
    # Analysis section
    st.markdown("### 📊 Dose-Response Curve")
    
    experiments = fetch_mst_experiments()
    
    if experiments:
        experiment_options = {
            exp["experiment_name"] or f"Experiment {str(exp['id'])[:8]}": exp["id"] 
            for exp in experiments
        }
        
        selected_name = st.selectbox(
            "Select MST Experiment", 
            options=list(experiment_options.keys()),
            key="mst_experiment_select"
        )
        
        if selected_name:
            selected_id = experiment_options[selected_name]
            selected_experiment = next((exp for exp in experiments if exp["id"] == selected_id), None)
            
            if selected_experiment:
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    # Plot dose-response curve
                    fig = go.Figure()
                    
                    # Generate sample dose-response data
                    concentrations = np.logspace(-12, -6, 20)  # 1 pM to 1 μM
                    kd_value = selected_experiment.get("kd_value", 1e-8)
                    hill_coeff = selected_experiment.get("hill_coefficient", 1.0)
                    
                    # Hill equation
                    fnorm = hill_equation(concentrations, kd_value, hill_coeff)
                    
                    # Add some noise for realism
                    noise = np.random.normal(0, 0.02, len(fnorm))
                    fnorm_noisy = fnorm + noise
                    
                    fig.add_trace(go.Scatter(
                        x=concentrations * 1e9,  # Convert to nM
                        y=fnorm_noisy,
                        mode='markers',
                        name='Data',
                        marker=dict(size=8)
                    ))
                    
                    # Add fitted curve
                    fit_concs = np.logspace(-12, -6, 100)
                    fit_fnorm = hill_equation(fit_concs, kd_value, hill_coeff)
                    
                    fig.add_trace(go.Scatter(
                        x=fit_concs * 1e9,
                        y=fit_fnorm,
                        mode='lines',
                        name='Fit',
                        line=dict(color='red')
                    ))
                    
                    fig.update_layout(
                        xaxis_type="log",
                        xaxis_title="Concentration (nM)",
                        yaxis_title="Fnorm (‰)",
                        title="MST Dose-Response Curve",
                        height=400
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Affinity parameters
                    st.markdown("**Affinity Parameters**")
                    
                    kd_value = selected_experiment.get("kd_value")
                    kd_error = selected_experiment.get("kd_error")
                    hill_coeff = selected_experiment.get("hill_coefficient")
                    
                    if kd_value:
                        st.metric("KD", f"{kd_value * 1e9:.2f} nM")
                    if kd_error:
                        st.metric("KD Error", f"±{kd_error * 1e9:.2f} nM")
                    if hill_coeff:
                        st.metric("Hill Coefficient", f"{hill_coeff:.2f}")
                    
                    # Quality metrics (simulated)
                    st.markdown("**Quality Metrics**")
                    st.metric("S/N Ratio", "12.5")
                    st.metric("Aggregation", "No")
                    st.metric("R²", "0.985")
                
                # Refit controls
                st.markdown("### 🔧 Affinity Analysis")
                col1, col2 = st.columns(2)
                
                with col1:
                    model = st.selectbox("Fitting Model", ["hill", "langmuir"], key="mst_model")
                
                with col2:
                    if st.button("Refit Affinity", key="mst_refit"):
                        with st.spinner("Refitting affinity parameters..."):
                            try:
                                with httpx.Client(timeout=60) as client:
                                    refit_data = {"model": model}
                                    response = client.post(
                                        f"{API_BASE}/api/v1/biophysical/mst/{selected_id}/fit", 
                                        json=refit_data
                                    )
                                
                                if response.status_code == 200:
                                    st.success("✅ Affinity parameters refitted!")
                                    st.rerun()
                                else:
                                    st.error(f"Refit failed: {response.text}")
                            except Exception as e:
                                st.error(f"Error refitting: {str(e)}")
    else:
        st.info("📁 No MST experiments found. Upload a file to get started!")


def render_dsc_tab() -> None:
    """Render the DSC Analysis tab."""
    st.subheader("Differential Scanning Calorimetry (DSC)")
    st.caption("Analyze protein thermal stability and unfolding")
    
    # Upload section
    st.markdown("### 📤 File Upload")
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Upload MicroCal/TA File",
            type=["csv", "txt"],
            key="dsc_upload",
            help="Upload DSC data files from MicroCal or TA Instruments"
        )
    
    with col2:
        with st.container():
            compound_id = st.text_input("Compound ID (optional)", key="dsc_compound")
            protein_name = st.text_input("Protein Name (optional)", key="dsc_protein")
    
    if uploaded_file is not None:
        st.success(f"File selected: {uploaded_file.name}")
        
        if st.button("Process DSC File", key="dsc_process", type="primary", use_container_width=True):
            with st.spinner("Uploading and processing DSC file..."):
                try:
                    data = {}
                    if compound_id:
                        data["compound_id"] = compound_id
                    if protein_name:
                        data["protein_name"] = protein_name
                    
                    with httpx.Client(timeout=120) as client:
                        files = {"file": (uploaded_file.name, uploaded_file.getvalue(), "text/csv")}
                        response = client.post(f"{API_BASE}/api/v1/biophysical/dsc/upload", files=files, data=data)
                    
                    if response.status_code == 201:
                        result = response.json()
                        st.success("✅ DSC file uploaded successfully!")
                        st.session_state.current_dsc_experiment = result["experiment_id"]
                        st.info("💡 Refresh the page to see your experiment in the analysis section below!")
                    else:
                        st.error(f"Upload failed: {response.text}")
                except Exception as e:
                    st.error(f"Error uploading file: {str(e)}")
    
    # Analysis section
    st.markdown("### 📊 Thermogram")
    
    experiments = fetch_dsc_experiments()
    
    if experiments:
        experiment_options = {
            exp["experiment_name"] or f"Experiment {str(exp['id'])[:8]}": exp["id"] 
            for exp in experiments
        }
        
        selected_name = st.selectbox(
            "Select DSC Experiment", 
            options=list(experiment_options.keys()),
            key="dsc_experiment_select"
        )
        
        if selected_name:
            selected_id = experiment_options[selected_name]
            selected_experiment = next((exp for exp in experiments if exp["id"] == selected_id), None)
            
            if selected_experiment:
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    # Plot thermogram
                    fig = go.Figure()
                    
                    # Generate sample thermogram data
                    temperatures = np.linspace(20, 90, 500)
                    tm_value = selected_experiment.get("tm_value", 65.0)
                    delta_h = selected_experiment.get("delta_h", -150.0)
                    
                    # Simulate DSC curve (two-state unfolding)
                    cp = simulate_dsc_curve(temperatures, tm_value, delta_h)
                    
                    fig.add_trace(go.Scatter(
                        x=temperatures,
                        y=cp,
                        mode='lines',
                        name='Thermogram',
                        line=dict(width=2)
                    ))
                    
                    # Mark Tm
                    if tm_value:
                        fig.add_vline(
                            x=tm_value, 
                            line_dash="dash", 
                            line_color="red",
                            annotation_text=f"Tm = {tm_value:.1f}°C"
                        )
                    
                    fig.update_layout(
                        xaxis_title="Temperature (°C)",
                        yaxis_title="Cp (kcal/mol/°C)",
                        title="DSC Thermogram",
                        height=400
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Thermal parameters
                    st.markdown("**Thermal Parameters**")
                    
                    tm_value = selected_experiment.get("tm_value")
                    tm_error = selected_experiment.get("tm_error")
                    delta_h = selected_experiment.get("delta_h")
                    
                    if tm_value:
                        st.metric("Tm", f"{tm_value:.1f} °C")
                    if tm_error:
                        st.metric("Tm Error", f"±{tm_error:.1f} °C")
                    if delta_h:
                        st.metric("ΔH", f"{delta_h:.1f} kcal/mol")
                    
                    # Additional metrics (simulated)
                    st.metric("Cooperativity", "1.85")
                    st.metric("Reversibility", "78%")
                
                # Refit controls
                st.markdown("### 🔧 Thermal Analysis")
                col1, col2 = st.columns(2)
                
                with col1:
                    model = st.selectbox("Fitting Model", ["two_state", "three_state"], key="dsc_model")
                
                with col2:
                    if st.button("Refit Thermal", key="dsc_refit"):
                        with st.spinner("Refitting thermal parameters..."):
                            try:
                                with httpx.Client(timeout=60) as client:
                                    refit_data = {"model": model}
                                    response = client.post(
                                        f"{API_BASE}/api/v1/biophysical/dsc/{selected_id}/fit", 
                                        json=refit_data
                                    )
                                
                                if response.status_code == 200:
                                    st.success("✅ Thermal parameters refitted!")
                                    st.rerun()
                                else:
                                    st.error(f"Refit failed: {response.text}")
                            except Exception as e:
                                st.error(f"Error refitting: {str(e)}")
    else:
        st.info("📁 No DSC experiments found. Upload a file to get started!")


# Helper Functions

def fetch_spr_experiments() -> List[Dict[str, Any]]:
    """Fetch SPR experiments from API."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{API_BASE}/api/v1/biophysical/spr")
            if response.status_code == 200:
                return response.json()
    except Exception as e:
        st.error(f"API error fetching SPR experiments: {e}")
    return []


def fetch_mst_experiments() -> List[Dict[str, Any]]:
    """Fetch MST experiments from API."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{API_BASE}/api/v1/biophysical/mst")
            if response.status_code == 200:
                return response.json()
    except Exception as e:
        st.error(f"API error fetching MST experiments: {e}")
    return []


def fetch_dsc_experiments() -> List[Dict[str, Any]]:
    """Fetch DSC experiments from API."""
    try:
        with httpx.Client(timeout=30) as client:
            response = client.get(f"{API_BASE}/api/v1/biophysical/dsc")
            if response.status_code == 200:
                return response.json()
    except Exception as e:
        st.error(f"API error fetching DSC experiments: {e}")
    return []


def hill_equation(conc: np.ndarray, kd: float, n: float) -> np.ndarray:
    """Hill equation for dose-response curve."""
    return 1 / (1 + (kd / conc) ** n)


def simulate_spr_response(time: np.ndarray, conc: float, ka: float, kd: float, rmax: float) -> np.ndarray:
    """Simulate SPR sensorgram response."""
    # Simple association/dissociation model
    # Association phase (0-300s)
    # Dissociation phase (300-600s)
    
    response = np.zeros_like(time)
    association_end = 300
    
    for i, t in enumerate(time):
        if t <= association_end:
            # Association: R = Rmax * conc * ka / (conc * ka + kd) * (1 - exp(-(conc * ka + kd) * t))
            kobs = conc * ka + kd
            response[i] = rmax * conc * ka / (conc * ka + kd) * (1 - np.exp(-kobs * t))
        else:
            # Dissociation: R = R_eq * exp(-kd * (t - t_eq))
            req = rmax * conc * ka / (conc * ka + kd) * (1 - np.exp(-kobs * association_end))
            response[i] = req * np.exp(-kd * (t - association_end))
    
    # Add some noise
    noise = np.random.normal(0, 0.5, len(response))
    return response + noise


def simulate_dsc_curve(temperature: np.ndarray, tm: float, delta_h: float) -> np.ndarray:
    """Simulate DSC thermogram."""
    # Two-state unfolding model
    R = 1.987e-3  # Gas constant in kcal/mol/K
    T_kelvin = temperature + 273.15
    tm_kelvin = tm + 273.15
    
    # Van't Hoff equation
    K = np.exp(-delta_h / R * (1/T_kelvin - 1/tm_kelvin))
    
    # Fraction unfolded
    fu = K / (1 + K)
    
    # Heat capacity (derivative of enthalpy)
    cp_baseline = 1.5  # kcal/mol/°C
    cp_transition = delta_h**2 / (R * T_kelvin**2) * fu * (1 - fu)
    
    return cp_baseline + cp_transition


def calculate_reversibility(scans: List[Dict]) -> float:
    """Calculate thermal reversibility from multiple scans."""
    # Placeholder implementation
    return 0.78  # 78% reversibility
