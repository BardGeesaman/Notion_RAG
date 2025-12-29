"""Experiment Planner page."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_experiment_planner_page() -> None:
    """Render the Experiment Planner page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📐 Experiment Planner")
    st.caption("Statistical power analysis, plate layouts, and cost estimation.")
    
    tab1, tab2, tab3 = st.tabs(["Power Calculator", "Plate Layout", "Cost Estimator"])
    
    with tab1:
        render_power_calculator_tab()
    
    with tab2:
        render_plate_layout_tab()
    
    with tab3:
        render_cost_estimator_tab()


def render_power_calculator_tab() -> None:
    """Render power calculator tab."""
    st.subheader("Statistical Power Calculator")
    
    # Effect size input
    use_preset = st.checkbox("Use preset effect size", value=True)
    
    if use_preset:
        effect_preset = st.selectbox(
            "Effect Size",
            options=["Small (0.2)", "Medium (0.5)", "Large (0.8)"],
            index=1,
        )
        effect_size = {"Small (0.2)": 0.2, "Medium (0.5)": 0.5, "Large (0.8)": 0.8}[effect_preset]
    else:
        effect_size = st.number_input("Custom Effect Size", min_value=0.01, max_value=5.0, value=0.5, step=0.1)
    
    col1, col2 = st.columns(2)
    
    with col1:
        alpha = st.number_input("Significance Level (α)", min_value=0.001, max_value=0.1, value=0.05, step=0.01)
    
    with col2:
        power = st.number_input("Desired Power", min_value=0.5, max_value=0.99, value=0.80, step=0.05)
    
    test_type = st.selectbox(
        "Test Type",
        options=["t-test", "anova", "correlation", "chi-square"],
        index=0,
    )
    
    if st.button("Calculate Sample Size", type="primary"):
        try:
            result = _api_post(
                "/api/v1/planner/power",
                {
                    "effect_size": effect_size,
                    "alpha": alpha,
                    "power": power,
                    "test_type": test_type,
                },
            )
            
            st.session_state["calculated_n"] = result.get("n_per_group")
            
        except Exception as e:
            st.error(f"Calculation failed: {e}")
    
    # Display result
    if "calculated_n" in st.session_state:
        n = st.session_state["calculated_n"]
        st.success(f"Required sample size: **{n} per group**")
        st.info(f"Total samples needed: {n * 2} (for 2 groups)")


def render_plate_layout_tab() -> None:
    """Render plate layout tab."""
    st.subheader("Plate Layout Calculator")
    
    # Sample size input
    n = st.number_input(
        "Number of Samples",
        min_value=1,
        max_value=10000,
        value=st.session_state.get("calculated_n", 96),
        step=1,
        help="Total samples to plate",
    )
    
    plate_format = st.selectbox(
        "Plate Format",
        options=[96, 384, 1536],
        index=0,
        help="Wells per plate",
    )
    
    if st.button("Calculate Layout", type="primary"):
        try:
            result = _api_post(
                "/api/v1/planner/plates",
                {"n": n, "plate_format": plate_format},
            )
            
            st.session_state["plate_layout"] = result
            
        except Exception as e:
            st.error(f"Calculation failed: {e}")
    
    # Display result
    layout = st.session_state.get("plate_layout")
    
    if layout:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Plates Needed", layout.get("plates_needed"))
        
        with col2:
            st.metric("Wells Used", layout.get("wells_used"))
        
        with col3:
            st.metric("Empty Wells", layout.get("empty_wells"))


def render_cost_estimator_tab() -> None:
    """Render cost estimator tab."""
    st.subheader("Experiment Cost Estimator")
    
    # Use calculated n if available
    default_n = st.session_state.get("calculated_n", 50)
    
    n = st.number_input(
        "Number of Samples",
        min_value=1,
        max_value=10000,
        value=default_n,
        step=1,
    )
    
    col1, col2 = st.columns(2)
    
    with col1:
        cost_per_sample = st.number_input(
            "Cost per Sample ($)",
            min_value=0.1,
            max_value=10000.0,
            value=10.0,
            step=1.0,
        )
    
    with col2:
        overhead_pct = st.slider(
            "Overhead %",
            min_value=0.0,
            max_value=50.0,
            value=10.0,
            step=5.0,
        ) / 100.0
    
    if st.button("Estimate Cost", type="primary"):
        try:
            result = _api_post(
                "/api/v1/planner/cost",
                {
                    "n": n,
                    "cost_per_sample": cost_per_sample,
                    "overhead_pct": overhead_pct,
                },
            )
            
            st.session_state["cost_estimate"] = result
            
        except Exception as e:
            st.error(f"Calculation failed: {e}")
    
    # Display result
    cost = st.session_state.get("cost_estimate")
    
    if cost:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Sample Cost", f"${cost.get('sample_cost', 0.0):,.2f}")
        
        with col2:
            st.metric("Overhead", f"${cost.get('overhead', 0.0):,.2f}")
        
        with col3:
            st.metric("Total Cost", f"${cost.get('total', 0.0):,.2f}")


if __name__ == "__main__":
    render_experiment_planner_page()

