"""
Dashboard widgets for Scientist's Cockpit.

Reusable widget components for overview dashboard.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 10) -> Any:
    """Make GET request to API."""
    try:
        with httpx.Client(timeout=timeout) as client:
            r = client.get(f"{API_BASE}{path}")
        r.raise_for_status()
        return r.json()
    except Exception as e:
        return None


def render_stats_widget() -> None:
    """
    Render statistics widget with entity counts.

    Shows counts for datasets, experiments, compounds, and signatures.
    """
    st.markdown("### 📊 Statistics")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Datasets count
        try:
            datasets = _api_get("/api/v1/datasets?limit=1")
            dataset_count = len(datasets) if datasets else 0
        except Exception:
            dataset_count = "?"
        
        st.metric("Datasets", dataset_count)
        
        # Compounds count
        compound_count = "?"  # Would need compounds endpoint
        st.metric("Compounds", compound_count)
    
    with col2:
        # Experiments count
        try:
            experiments = _api_get("/api/v1/experiments?limit=1")
            exp_count = len(experiments) if experiments else 0
        except Exception:
            exp_count = "?"
        
        st.metric("Experiments", exp_count)
        
        # Signatures count
        signature_count = "?"  # Would need signatures endpoint
        st.metric("Signatures", signature_count)
    
    if st.button("View All Data", key="stats_view_all", use_container_width=True):
        st.switch_page("pages/datasets.py")


def render_activity_widget() -> None:
    """
    Render recent activity widget.

    Shows last 5 activity events with timestamps.
    """
    st.markdown("### 📋 Recent Activity")
    
    try:
        events = _api_get("/api/v1/activity/feed?limit=5")
        
        if events:
            for event in events[:5]:
                event_type = event.get("event_type", "unknown")
                target_name = event.get("target_name", "Unknown")
                created_at = event.get("created_at", "")
                
                # Simple event display
                icon = "🔬" if "experiment" in event_type else "📊"
                st.markdown(f"{icon} {target_name}")
                if created_at:
                    st.caption(created_at[:10])  # Show date
        else:
            st.info("No recent activity")
    
    except Exception:
        st.warning("Could not load activity")
    
    if st.button("View Activity Feed", key="activity_view_all", use_container_width=True):
        st.switch_page("pages/activity_feed.py")


def render_alerts_widget() -> None:
    """
    Render alerts widget.

    Shows unread notifications with severity badges.
    """
    st.markdown("### 🔔 Alerts")
    
    try:
        # Try to get notification count
        notif_data = _api_get("/api/v1/notifications/count")
        
        if notif_data:
            unread = notif_data.get("unread_count", 0)
            
            if unread > 0:
                st.warning(f"{unread} unread notifications")
            else:
                st.success("No new alerts")
        else:
            st.info("No alerts")
    
    except Exception:
        st.info("Alerts unavailable")
    
    if st.button("View All Alerts", key="alerts_view_all", use_container_width=True):
        st.switch_page("pages/overview.py")


def render_tasks_widget() -> None:
    """
    Render pending tasks widget.

    Shows pending reviews and approvals.
    """
    st.markdown("### ✅ Tasks")
    
    try:
        # Try to get pending reviews
        reviews = _api_get("/api/reviews/pending")
        
        if reviews:
            pending_count = len(reviews) if isinstance(reviews, list) else 0
            
            if pending_count > 0:
                st.warning(f"{pending_count} pending reviews")
                
                # Show first few
                for review in reviews[:3]:
                    notebook_path = review.get("notebook_path", "Unknown")
                    st.caption(f"📝 {notebook_path}")
            else:
                st.success("No pending tasks")
        else:
            st.info("No pending tasks")
    
    except Exception:
        st.info("Tasks unavailable")
    
    if st.button("View Review Queue", key="tasks_view_all", use_container_width=True):
        st.switch_page("pages/review_queue.py")


def render_shortcuts_widget() -> None:
    """
    Render quick action shortcuts widget.

    Provides buttons for common actions.
    """
    st.markdown("### ⚡ Quick Actions")
    
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("➕ New Experiment", key="shortcut_experiment", use_container_width=True):
            st.switch_page("pages/experiments.py")
        
        if st.button("📥 Upload Data", key="shortcut_upload", use_container_width=True):
            st.switch_page("pages/import_data.py")
    
    with col2:
        if st.button("🔍 Search Papers", key="shortcut_papers", use_container_width=True):
            st.switch_page("pages/paper_search.py")
        
        if st.button("🧪 HTS QC", key="shortcut_hts", use_container_width=True):
            st.switch_page("pages/hts_qc.py")

