"""Share Links Dashboard page."""

from __future__ import annotations

import os
from typing import Any

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


def _api_post(path: str, payload: dict, *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_delete(path: str, *, timeout: int = 30) -> Any:
    """Make DELETE request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.delete(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def render_share_links_page() -> None:
    """Render the Share Links page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🔗 Share Links")
    st.caption("Create shareable links for Voila dashboards.")
    
    tab1, tab2 = st.tabs(["Create Share Link", "Manage Links"])
    
    with tab1:
        render_create_tab()
    
    with tab2:
        render_manage_tab()


def render_create_tab() -> None:
    """Render create share link tab."""
    st.subheader("Generate Shareable Link")
    
    dashboard_path = st.text_input(
        "Dashboard Path",
        placeholder="/voila/render/templates/hts_plate_viewer.ipynb",
        help="Path to Voila dashboard notebook",
    )
    
    context_json = st.text_area(
        "Context JSON (optional)",
        placeholder='{"experiment_id": "123", "campaign_id": "456"}',
        height=100,
        help="Dashboard context parameters as JSON",
    )
    
    col1, col2 = st.columns(2)
    
    with col1:
        expiration = st.selectbox(
            "Expires In",
            options=["1 hour", "24 hours", "7 days", "30 days"],
            index=1,
        )
        
        # Convert to hours
        expires_map = {
            "1 hour": 1,
            "24 hours": 24,
            "7 days": 168,
            "30 days": 720,
        }
        expires_in_hours = expires_map[expiration]
    
    with col2:
        max_views = st.number_input(
            "Max Views (0 = unlimited)",
            min_value=0,
            max_value=1000,
            value=0,
            step=1,
        )
    
    permissions = st.radio(
        "Permissions",
        options=["View Only", "Interactive"],
        horizontal=True,
    )
    
    # User ID for auth (would come from session in production)
    user_id = st.text_input("User ID", placeholder="Your user UUID")
    company_id = st.text_input("Company ID", placeholder="Your company UUID")
    
    if st.button("🔗 Generate Link", type="primary", disabled=not all([dashboard_path, user_id, company_id])):
        try:
            # Parse context if provided
            import json
            context = json.loads(context_json) if context_json else None
            
            result = _api_post(
                "/api/v1/share-links",
                {
                    "dashboard_path": dashboard_path,
                    "context": context,
                    "expires_in_hours": expires_in_hours,
                    "max_views": max_views if max_views > 0 else None,
                    "permissions": "view" if permissions == "View Only" else "interactive",
                },
            )
            
            st.success("✅ Share link created!")
            
            # Display shareable URL
            share_url = f"{API_BASE}/share/{result.get('token')}"
            st.code(share_url, language=None)
            
            st.download_button(
                "📋 Copy to Clipboard",
                data=share_url,
                file_name="share_link.txt",
            )
            
        except json.JSONDecodeError:
            st.error("Invalid JSON in context field")
        except Exception as e:
            st.error(f"Failed to create share link: {e}")


def render_manage_tab() -> None:
    """Render manage links tab."""
    st.subheader("Manage Share Links")
    
    if st.button("🔄 Load My Links"):
        try:
            links = _api_get("/api/v1/share-links")
            st.session_state["share_links"] = links
        except Exception as e:
            st.error(f"Failed to load links: {e}")
    
    links = st.session_state.get("share_links")
    
    if links:
        for link in links:
            with st.expander(f"📊 {link.get('dashboard_path', 'Unknown')}"):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Views", f"{link.get('view_count', 0)}/{link.get('max_views', '∞')}")
                
                with col2:
                    status = "🟢 Active" if link.get("is_active") else "🔴 Inactive"
                    st.metric("Status", status)
                
                with col3:
                    st.metric("Expires", link.get("expires_at", "")[:10])
                
                st.code(f"Token: {link.get('token', '')[:32]}...", language=None)
                
                col_a, col_b = st.columns(2)
                
                with col_a:
                    if st.button("📊 Stats", key=f"stats_{link.get('id')}"):
                        try:
                            stats = _api_get(f"/api/v1/share-links/{link.get('id')}/stats")
                            st.json(stats)
                        except Exception as e:
                            st.error(f"Failed to get stats: {e}")
                
                with col_b:
                    if st.button("🗑️ Revoke", key=f"revoke_{link.get('id')}"):
                        try:
                            _api_delete(f"/api/v1/share-links/{link.get('id')}")
                            st.success("Link revoked")
                            st.session_state.pop("share_links", None)
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to revoke: {e}")
    else:
        st.info("No share links yet. Create one in the Create Share Link tab!")


if __name__ == "__main__":
    render_share_links_page()

