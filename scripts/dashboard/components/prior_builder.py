"""Prior Builder widget for Bayesian dose-response modeling."""

from typing import Dict, Optional

import streamlit as st


def render_prior_builder(key_prefix: str = "prior") -> Optional[Dict]:
    """Render Prior Builder widget. Returns prior config dict or None.
    
    Args:
        key_prefix: Unique prefix for widget keys to avoid conflicts
        
    Returns:
        Dictionary with prior configuration or None
    """
    
    st.subheader("Prior Configuration")
    
    preset = st.selectbox(
        "Load Preset", 
        ["Custom", "Biochemical Assay", "Cell Viability", "Weak Prior"],
        key=f"{key_prefix}_preset"
    )
    
    # Preset defaults
    presets = {
        "Custom": {"ec50": 0.0, "ec50_sd": 1.0, "hill": 1.0, "hill_sd": 2.0},
        "Biochemical Assay": {"ec50": -6.0, "ec50_sd": 1.5, "hill": 1.0, "hill_sd": 1.0},
        "Cell Viability": {"ec50": -5.0, "ec50_sd": 2.0, "hill": 1.5, "hill_sd": 1.5},
        "Weak Prior": {"ec50": 0.0, "ec50_sd": 3.0, "hill": 1.0, "hill_sd": 3.0},
    }
    defaults = presets[preset]
    
    col1, col2 = st.columns(2)
    with col1:
        ec50_mean = st.slider(
            "EC50 Mean (log)", 
            -5.0, 5.0, 
            defaults["ec50"], 
            0.1, 
            key=f"{key_prefix}_ec50"
        )
        ec50_sd = st.slider(
            "EC50 SD", 
            0.1, 3.0, 
            defaults["ec50_sd"], 
            0.1, 
            key=f"{key_prefix}_ec50_sd"
        )
    with col2:
        hill_mean = st.slider(
            "Hill Mean", 
            -2.0, 4.0, 
            defaults["hill"], 
            0.1, 
            key=f"{key_prefix}_hill"
        )
        hill_sd = st.slider(
            "Hill SD", 
            0.1, 3.0, 
            defaults["hill_sd"], 
            0.1, 
            key=f"{key_prefix}_hill_sd"
        )
    
    return {
        "ec50_prior_mean": ec50_mean,
        "ec50_prior_sd": ec50_sd,
        "hill_prior_mean": hill_mean,
        "hill_prior_sd": hill_sd,
    }
