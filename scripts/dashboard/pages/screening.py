"""HTS Screening Campaigns page."""

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


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_screening_page() -> None:
    """Render the HTS Screening Campaigns page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🧪 HTS Screening Campaigns")
    st.caption("View high-throughput screening campaigns, hits, and active learning suggestions.")
    
    tab1, tab2, tab3 = st.tabs(["Campaigns", "Campaign Details", "Active Learning"])
    
    with tab1:
        render_campaigns_tab()
    
    with tab2:
        render_campaign_details_tab()
    
    with tab3:
        render_active_learning_tab()


def render_campaigns_tab() -> None:
    """Render the campaigns list tab."""
    st.subheader("Screening Campaigns")
    
    try:
        campaigns = _api_get("/api/screening/campaigns")
    except Exception as e:
        st.error(f"Failed to load campaigns: {e}")
        return
    
    if not campaigns:
        st.info("No screening campaigns found.")
        return
    
    # Display campaigns table
    campaign_data = []
    for c in campaigns:
        campaign_data.append({
            "Campaign ID": c.get("campaign_id", ""),
            "Name": c.get("campaign_name", ""),
            "Assay Type": c.get("assay_type", ""),
            "Target": c.get("target", ""),
            "Total Wells": c.get("total_wells", 0),
            "Hits": c.get("hit_count", 0),
            "Hit Rate %": f"{(c.get('hit_count', 0) / c.get('total_wells', 1) * 100):.1f}" if c.get('total_wells', 0) > 0 else "0.0",
            "Run Date": c.get("run_date", ""),
        })
    
    df = pd.DataFrame(campaign_data)
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    # Store selected campaign
    if len(campaigns) > 0:
        campaign_ids = [c.get("campaign_id", "") for c in campaigns if c.get("campaign_id")]
        if campaign_ids:
            selected = st.selectbox("Select campaign for details", [""] + campaign_ids)
            if selected:
                st.session_state["selected_campaign_id"] = selected


def render_campaign_details_tab() -> None:
    """Render the campaign details tab."""
    st.subheader("Campaign Details")
    
    campaign_id = st.session_state.get("selected_campaign_id")
    
    if not campaign_id:
        campaign_id = st.text_input("Campaign ID", placeholder="Enter campaign ID...")
    
    if not campaign_id:
        st.info("Select a campaign from the Campaigns tab or enter a campaign ID above.")
        return
    
    # Fetch campaign details
    try:
        campaign = _api_get(f"/api/screening/campaigns/{campaign_id}")
    except httpx.HTTPStatusError as e:
        if e.response.status_code == 404:
            st.error(f"Campaign '{campaign_id}' not found.")
        else:
            st.error(f"Failed to load campaign: {e}")
        return
    except Exception as e:
        st.error(f"Failed to load campaign: {e}")
        return
    
    # Display campaign metadata
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Campaign ID", campaign.get("campaign_id", ""))
    with col2:
        st.metric("Total Wells", campaign.get("total_wells", 0))
    with col3:
        st.metric("Hits", campaign.get("hit_count", 0))
    
    st.markdown(f"**Name**: {campaign.get('campaign_name', 'N/A')}")
    st.markdown(f"**Target**: {campaign.get('target', 'N/A')}")
    st.markdown(f"**Assay Type**: {campaign.get('assay_type', 'N/A')}")
    
    if campaign.get("description"):
        st.markdown(f"**Description**: {campaign.get('description')}")
    
    # Fetch and display hits
    st.markdown("---")
    st.subheader("Campaign Hits")
    
    try:
        hits = _api_get(f"/api/screening/campaigns/{campaign_id}/hits")
    except Exception as e:
        st.error(f"Failed to load hits: {e}")
        return
    
    if not hits:
        st.info("No hits found for this campaign.")
        return
    
    # Display hits table
    hits_data = []
    for h in hits:
        hits_data.append({
            "Compound ID": h.get("compound_id", ""),
            "Well": h.get("well_position", ""),
            "Raw Value": h.get("raw_value", ""),
            "Normalized": h.get("normalized_value", ""),
            "Z-Score": h.get("z_score", ""),
            "Category": h.get("hit_category", ""),
        })
    
    df_hits = pd.DataFrame(hits_data)
    st.dataframe(df_hits, use_container_width=True, hide_index=True)
    
    # Download hits
    csv = df_hits.to_csv(index=False)
    st.download_button(
        "📥 Download Hits (CSV)",
        data=csv,
        file_name=f"hits_{campaign_id}.csv",
        mime="text/csv",
    )


def render_active_learning_tab() -> None:
    """Render the active learning suggestions tab."""
    st.subheader("Active Learning Suggestions")
    st.caption("Get AI-powered compound suggestions for next screening round.")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        strategy = st.selectbox(
            "Strategy",
            options=["uncertainty", "diversity", "exploitation"],
            index=0,
            help="Acquisition function for active learning",
        )
    
    with col2:
        batch_size = st.number_input(
            "Batch Size",
            min_value=1,
            max_value=100,
            value=10,
            step=1,
            help="Number of compounds to suggest",
        )
    
    model_id = st.text_input(
        "Model ID (optional)",
        placeholder="UUID of trained model for uncertainty sampling...",
        help="Required for uncertainty strategy",
    )
    
    st.markdown("**Note**: This is a placeholder UI. Full active learning integration requires:")
    st.markdown("- Uploading screened compounds")
    st.markdown("- Uploading candidate pool")
    st.markdown("- Training a model (for uncertainty sampling)")
    
    if st.button("Get Suggestions", type="primary"):
        st.info("Active learning suggestion API ready. Upload screened compounds and candidates to generate suggestions.")


if __name__ == "__main__":
    render_screening_page()

