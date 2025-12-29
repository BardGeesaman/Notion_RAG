"""Compound Portfolio Dashboard page."""

from __future__ import annotations

import os
from typing import Any

import httpx
import plotly.express as px
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def render_compound_portfolio_page() -> None:
    """Render the Compound Portfolio Dashboard page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("💼 Compound Portfolio")
    st.caption("Unified view of compound collection with ADMET, SAR, and recommendations.")
    
    tab1, tab2, tab3, tab4 = st.tabs(["Overview", "ADMET Summary", "SAR Gaps", "Recommendations"])
    
    with tab1:
        render_overview_tab()
    
    with tab2:
        render_admet_tab()
    
    with tab3:
        render_sar_gaps_tab()
    
    with tab4:
        render_recommendations_tab()


def render_overview_tab() -> None:
    """Render portfolio overview tab."""
    st.subheader("Portfolio Overview")
    
    try:
        summary = _api_get("/api/v1/portfolio/summary")
        
        # Metrics cards
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Compounds", summary.get("total_compounds", 0))
        
        with col2:
            st.metric("Scaffolds", summary.get("scaffold_count", 0))
        
        with col3:
            date_from = summary.get("date_from", "")
            st.metric("First Compound", date_from[:10] if date_from else "N/A")
        
        with col4:
            date_to = summary.get("date_to", "")
            st.metric("Latest Compound", date_to[:10] if date_to else "N/A")
    
    except Exception as e:
        st.error(f"Failed to load portfolio summary: {e}")


def render_admet_tab() -> None:
    """Render ADMET summary tab."""
    st.subheader("ADMET Profile")
    
    try:
        admet = _api_get("/api/v1/portfolio/admet")
        
        # Traffic light counts
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("✅ Green", admet.get("green", 0), help="ADMET-compliant")
        
        with col2:
            st.metric("⚠️ Yellow", admet.get("yellow", 0), help="Borderline")
        
        with col3:
            st.metric("🚫 Red", admet.get("red", 0), help="Failed ADMET")
        
        with col4:
            st.metric("❓ Unknown", admet.get("unknown", 0), help="Not predicted")
        
        # Pie chart
        if any([admet.get("green"), admet.get("yellow"), admet.get("red")]):
            import pandas as pd
            
            df = pd.DataFrame({
                "Status": ["Green", "Yellow", "Red", "Unknown"],
                "Count": [
                    admet.get("green", 0),
                    admet.get("yellow", 0),
                    admet.get("red", 0),
                    admet.get("unknown", 0),
                ],
            })
            
            fig = px.pie(
                df,
                values="Count",
                names="Status",
                title="ADMET Distribution",
                color="Status",
                color_discrete_map={
                    "Green": "green",
                    "Yellow": "yellow",
                    "Red": "red",
                    "Unknown": "gray",
                },
            )
            st.plotly_chart(fig, use_container_width=True)
    
    except Exception as e:
        st.error(f"Failed to load ADMET summary: {e}")


def render_sar_gaps_tab() -> None:
    """Render SAR gaps tab."""
    st.subheader("SAR Coverage Gaps")
    st.caption("Scaffolds with sparse compound coverage")
    
    try:
        gaps = _api_get("/api/v1/portfolio/gaps")
        
        if gaps:
            import pandas as pd
            
            df = pd.DataFrame(gaps)
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            st.info(f"Found {len(gaps)} scaffolds needing expansion (< 3 compounds)")
        else:
            st.success("No significant SAR gaps detected!")
    
    except Exception as e:
        st.error(f"Failed to load SAR gaps: {e}")


def render_recommendations_tab() -> None:
    """Render recommendations tab."""
    st.subheader("Next Compound Suggestions")
    st.caption("Top-ranked compounds for synthesis or testing")
    
    try:
        recommendations = _api_get("/api/v1/portfolio/recommendations?limit=10")
        
        if recommendations:
            for i, rec in enumerate(recommendations, 1):
                with st.expander(f"{i}. {rec.get('compound_id', 'Unknown')} (Score: {rec.get('score', 0.0):.2f})"):
                    st.code(rec.get("smiles", ""), language=None)
                    st.markdown(f"**Reason:** {rec.get('reason', 'N/A')}")
            
            st.markdown("---")
            if st.button("➕ Draw New Compound"):
                st.switch_page("pages/chemical_sketcher.py")
        else:
            st.info("No recommendations available. Add more compounds to enable suggestions.")
    
    except Exception as e:
        st.error(f"Failed to load recommendations: {e}")


if __name__ == "__main__":
    render_compound_portfolio_page()

