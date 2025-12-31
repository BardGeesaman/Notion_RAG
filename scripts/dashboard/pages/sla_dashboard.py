"""SLA Dashboard - Monitor and manage SLA rules and review cycles."""

import os
from datetime import datetime, timezone
from typing import Dict

import httpx
import streamlit as st

from scripts.dashboard.core.auth import check_authentication
from scripts.dashboard.core.config import AUTH_DISABLED

API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str) -> Dict:
    """Make authenticated GET request to API."""
    try:
        response = httpx.get(f"{API_BASE}{path}", timeout=10.0)
        response.raise_for_status()
        return response.json()
    except httpx.ConnectError:
        raise
    except Exception as e:
        st.error(f"API request failed: {e}")
        raise


def _api_post(path: str, data: Dict) -> Dict:
    """Make authenticated POST request to API."""
    try:
        response = httpx.post(f"{API_BASE}{path}", json=data, timeout=10.0)
        response.raise_for_status()
        return response.json()
    except httpx.ConnectError:
        raise
    except Exception as e:
        st.error(f"API request failed: {e}")
        raise


def _api_patch(path: str, data: Dict) -> Dict:
    """Make authenticated PATCH request to API."""
    try:
        response = httpx.patch(f"{API_BASE}{path}", json=data, timeout=10.0)
        response.raise_for_status()
        return response.json()
    except httpx.ConnectError:
        raise
    except Exception as e:
        st.error(f"API request failed: {e}")
        raise


def _api_delete(path: str) -> None:
    """Make authenticated DELETE request to API."""
    try:
        response = httpx.delete(f"{API_BASE}{path}", timeout=10.0)
        response.raise_for_status()
    except httpx.ConnectError:
        raise
    except Exception as e:
        st.error(f"API request failed: {e}")
        raise


def _format_time_remaining(due_at_str: str) -> str:
    """Format time remaining until due date."""
    try:
        due_at = datetime.fromisoformat(due_at_str.replace("Z", "+00:00"))
        now = datetime.now(timezone.utc)
        diff = due_at - now
        
        if diff.total_seconds() < 0:
            hours_overdue = abs(diff.total_seconds()) / 3600
            return f"Overdue by {hours_overdue:.1f}h"
        else:
            hours_remaining = diff.total_seconds() / 3600
            if hours_remaining < 24:
                return f"{hours_remaining:.1f}h remaining"
            else:
                days_remaining = hours_remaining / 24
                return f"{days_remaining:.1f}d remaining"
    except Exception:
        return "Unknown"


def _get_status_color(status: str) -> str:
    """Get color for SLA status."""
    colors = {
        "on_track": "🟢",
        "warning": "🟡", 
        "overdue": "🟠",
        "breached": "🔴",
    }
    return colors.get(status, "⚪")


def render_overview_tab():
    """Render the Overview tab with SLA metrics."""
    st.subheader("📊 SLA Status Overview")
    
    try:
        with st.spinner("Loading SLA status..."):
            status = _api_get("/api/v1/sla/status")
        
        # Metrics cards
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                "🟢 On Track",
                status.get("on_track", 0),
                help="Reviews progressing within SLA"
            )
        
        with col2:
            st.metric(
                "🟡 At Risk",
                status.get("warning", 0),
                help="Reviews approaching SLA deadline"
            )
        
        with col3:
            st.metric(
                "🟠 Overdue",
                status.get("overdue", 0),
                help="Reviews past SLA deadline"
            )
        
        with col4:
            st.metric(
                "🔴 Breached",
                status.get("breached", 0),
                help="Reviews significantly past SLA"
            )
        
        # Summary
        total = status.get("total", 0)
        if total > 0:
            on_track_pct = (status.get("on_track", 0) / total) * 100
            st.info(f"**Total Active Reviews**: {total} | **On Track**: {on_track_pct:.1f}%")
        else:
            st.info("No active reviews with SLAs")
            
    except httpx.ConnectError:
        st.warning("⚠️ API unavailable. Showing demo data.")
        # Demo data
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("🟢 On Track", 15)
        with col2:
            st.metric("🟡 At Risk", 3)
        with col3:
            st.metric("🟠 Overdue", 2)
        with col4:
            st.metric("🔴 Breached", 1)
        st.info("**Total Active Reviews**: 21 | **On Track**: 71.4%")
    except Exception as e:
        st.error(f"Failed to load SLA status: {e}")


def render_reviews_tab():
    """Render the Reviews tab with active reviews."""
    st.subheader("📋 Active Reviews")
    
    # Filters
    col1, col2 = st.columns(2)
    with col1:
        entity_filter = st.selectbox(
            "Entity Type",
            ["All", "dataset", "experiment", "compound", "signature"],
            key="reviews_entity_filter"
        )
    with col2:
        st.selectbox(
            "Status",
            ["All", "on_track", "warning", "overdue", "breached"],
            key="reviews_status_filter"
        )
    
    try:
        with st.spinner("Loading reviews..."):
            reviews = _api_get("/api/v1/sla/overdue?limit=100")
        
        if not reviews:
            st.info("No overdue reviews found.")
            return
        
        # Filter reviews
        filtered_reviews = reviews
        if entity_filter != "All":
            filtered_reviews = [r for r in filtered_reviews if r.get("entity_type") == entity_filter]
        
        # Display reviews table
        st.write(f"Found {len(filtered_reviews)} reviews")
        
        for review in filtered_reviews:
            with st.expander(f"{_get_status_color(review.get('sla_status', 'unknown'))} {review.get('entity_type', 'Unknown')} - {review.get('entity_id', 'Unknown')}"):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.write(f"**Entity Type**: {review.get('entity_type', 'Unknown')}")
                    st.write(f"**Entity ID**: {review.get('entity_id', 'Unknown')}")
                
                with col2:
                    st.write(f"**Reviewer**: {review.get('reviewer_id', 'Unassigned')}")
                    st.write(f"**Status**: {review.get('status', 'Unknown')}")
                
                with col3:
                    st.write(f"**SLA Status**: {_get_status_color(review.get('sla_status', 'unknown'))} {review.get('sla_status', 'Unknown').title()}")
                    if review.get("due_at"):
                        st.write(f"**Time**: {_format_time_remaining(review['due_at'])}")
                
    except httpx.ConnectError:
        st.warning("⚠️ API unavailable. Showing demo data.")
        # Demo data
        demo_reviews = [
            {
                "entity_type": "dataset",
                "entity_id": "demo-dataset-1",
                "reviewer_id": "alice@company.com",
                "status": "in_review",
                "sla_status": "warning",
                "due_at": "2025-01-01T18:00:00Z"
            },
            {
                "entity_type": "experiment",
                "entity_id": "demo-experiment-1", 
                "reviewer_id": "bob@company.com",
                "status": "pending",
                "sla_status": "overdue",
                "due_at": "2024-12-30T12:00:00Z"
            }
        ]
        
        for review in demo_reviews:
            with st.expander(f"{_get_status_color(review['sla_status'])} {review['entity_type']} - {review['entity_id']}"):
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.write(f"**Entity Type**: {review['entity_type']}")
                    st.write(f"**Entity ID**: {review['entity_id']}")
                with col2:
                    st.write(f"**Reviewer**: {review['reviewer_id']}")
                    st.write(f"**Status**: {review['status']}")
                with col3:
                    st.write(f"**SLA Status**: {_get_status_color(review['sla_status'])} {review['sla_status'].title()}")
                    st.write(f"**Time**: {_format_time_remaining(review['due_at'])}")
    except Exception as e:
        st.error(f"Failed to load reviews: {e}")


def render_cycles_tab(user):
    """Render the Cycles tab for managing review cycles."""
    st.subheader("🔄 Review Cycles")
    
    if not user or user.get("role") != "admin":
        st.warning("🔒 Admin access required to manage review cycles.")
        return
    
    # Create cycle form
    with st.expander("➕ Create New Cycle"):
        with st.form("create_cycle_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                cycle_name = st.text_input("Cycle Name", placeholder="e.g., Weekly Dataset Reviews")
                entity_type = st.selectbox("Entity Type", ["dataset", "experiment", "compound", "signature"])
                frequency = st.selectbox("Frequency", ["weekly", "monthly", "quarterly", "yearly"])
            
            with col2:
                if frequency == "weekly":
                    day_of_week = st.selectbox("Day of Week", {
                        0: "Monday", 1: "Tuesday", 2: "Wednesday", 3: "Thursday",
                        4: "Friday", 5: "Saturday", 6: "Sunday"
                    })
                    day_of_month = None
                else:
                    day_of_week = None
                    day_of_month = st.number_input("Day of Month", min_value=1, max_value=28, value=1)
                
                reviewer_pool_text = st.text_area(
                    "Reviewer Pool (one email per line)",
                    placeholder="alice@company.com\nbob@company.com"
                )
            
            if st.form_submit_button("Create Cycle", type="primary"):
                if cycle_name and reviewer_pool_text:
                    reviewer_pool = [email.strip() for email in reviewer_pool_text.split("\n") if email.strip()]
                    
                    cycle_data = {
                        "name": cycle_name,
                        "entity_type": entity_type,
                        "frequency": frequency,
                        "reviewer_pool": reviewer_pool,
                    }
                    
                    if day_of_week is not None:
                        cycle_data["day_of_week"] = day_of_week
                    if day_of_month is not None:
                        cycle_data["day_of_month"] = day_of_month
                    
                    try:
                        with st.spinner("Creating cycle..."):
                            _api_post("/api/v1/sla/review-cycles", cycle_data)
                        st.success("✅ Cycle created successfully!")
                        st.rerun()
                    except httpx.ConnectError:
                        st.error("⚠️ API unavailable. Cannot create cycle.")
                    except Exception as e:
                        st.error(f"Failed to create cycle: {e}")
                else:
                    st.warning("Please fill in all required fields.")
    
    # List existing cycles
    st.divider()
    
    try:
        with st.spinner("Loading cycles..."):
            cycles = _api_get("/api/v1/sla/review-cycles")
        
        if not cycles:
            st.info("No review cycles configured.")
            return
        
        st.write(f"**{len(cycles)} Active Cycles**")
        
        for cycle in cycles:
            with st.expander(f"{'✅' if cycle.get('is_active') else '❌'} {cycle.get('name', 'Unnamed')} ({cycle.get('frequency', 'Unknown')})"):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.write(f"**Entity Type**: {cycle.get('entity_type', 'Unknown')}")
                    st.write(f"**Frequency**: {cycle.get('frequency', 'Unknown')}")
                    if cycle.get("day_of_week") is not None:
                        days = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
                        st.write(f"**Day of Week**: {days[cycle['day_of_week']]}")
                    if cycle.get("day_of_month"):
                        st.write(f"**Day of Month**: {cycle['day_of_month']}")
                
                with col2:
                    st.write(f"**Active**: {'Yes' if cycle.get('is_active') else 'No'}")
                    if cycle.get("next_run_at"):
                        next_run = datetime.fromisoformat(cycle["next_run_at"].replace("Z", "+00:00"))
                        st.write(f"**Next Run**: {next_run.strftime('%Y-%m-%d %H:%M UTC')}")
                
                with col3:
                    reviewer_count = len(cycle.get("reviewer_pool", []))
                    st.write(f"**Reviewers**: {reviewer_count}")
                    
                    # Action buttons
                    col_run, col_edit, col_delete = st.columns(3)
                    
                    with col_run:
                        if st.button("▶️ Run Now", key=f"run_{cycle['id']}"):
                            try:
                                with st.spinner("Running cycle..."):
                                    result = _api_post(f"/api/v1/sla/review-cycles/{cycle['id']}/run-now", {})
                                st.success(f"✅ {result.get('message', 'Cycle executed')}")
                            except Exception as e:
                                st.error(f"Failed to run cycle: {e}")
                    
                    with col_edit:
                        if st.button("✏️ Edit", key=f"edit_{cycle['id']}"):
                            st.info("Edit functionality coming soon")
                    
                    with col_delete:
                        if st.button("🗑️ Delete", key=f"delete_{cycle['id']}"):
                            try:
                                with st.spinner("Deleting cycle..."):
                                    _api_delete(f"/api/v1/sla/review-cycles/{cycle['id']}")
                                st.success("✅ Cycle deleted!")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Failed to delete cycle: {e}")
        
    except httpx.ConnectError:
        st.warning("⚠️ API unavailable. Showing demo data.")
        # Demo data
        st.write("**2 Active Cycles**")
        
        with st.expander("✅ Weekly Dataset Reviews (weekly)"):
            st.write("**Entity Type**: dataset")
            st.write("**Frequency**: weekly")
            st.write("**Next Run**: 2025-01-06 09:00 UTC")
            st.write("**Reviewers**: 3")
    except Exception as e:
        st.error(f"Failed to load cycles: {e}")


def render_settings_tab(user):
    """Render the Settings tab for managing SLA rules."""
    st.subheader("⚙️ SLA Rules")
    
    if not user or user.get("role") != "admin":
        st.warning("🔒 Admin access required to manage SLA rules.")
        return
    
    # Create SLA form
    with st.expander("➕ Create New SLA Rule"):
        with st.form("create_sla_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                sla_name = st.text_input("SLA Name", placeholder="e.g., Critical Dataset SLA")
                entity_type = st.selectbox("Entity Type", ["All Types", "dataset", "experiment", "compound", "signature"])
                max_hours = st.number_input("Max Review Hours", min_value=1, max_value=720, value=120)
            
            with col2:
                warning_threshold = st.slider("Warning Threshold (%)", min_value=50, max_value=95, value=75)
                is_default = st.checkbox("Set as Default SLA")
                escalation_emails = st.text_area(
                    "Escalation Chain (one email per line)",
                    placeholder="manager@company.com\ndirector@company.com"
                )
            
            if st.form_submit_button("Create SLA Rule", type="primary"):
                if sla_name:
                    escalation_chain = None
                    if escalation_emails:
                        escalation_chain = [email.strip() for email in escalation_emails.split("\n") if email.strip()]
                    
                    sla_data = {
                        "name": sla_name,
                        "entity_type": None if entity_type == "All Types" else entity_type,
                        "max_review_hours": max_hours,
                        "warning_threshold_pct": warning_threshold,
                        "escalation_chain": escalation_chain,
                        "is_default": is_default,
                        "is_active": True,
                    }
                    
                    try:
                        with st.spinner("Creating SLA rule..."):
                            _api_post("/api/v1/sla/rules", sla_data)
                        st.success("✅ SLA rule created successfully!")
                        st.rerun()
                    except httpx.ConnectError:
                        st.error("⚠️ API unavailable. Cannot create SLA rule.")
                    except Exception as e:
                        st.error(f"Failed to create SLA rule: {e}")
                else:
                    st.warning("Please enter an SLA name.")
    
    # List existing SLA rules
    st.divider()
    
    try:
        with st.spinner("Loading SLA rules..."):
            sla_rules = _api_get("/api/v1/sla/rules")
        
        if not sla_rules:
            st.info("No SLA rules configured.")
            return
        
        st.write(f"**{len(sla_rules)} SLA Rules**")
        
        for sla in sla_rules:
            entity_display = sla.get("entity_type", "All Types") or "All Types"
            default_badge = "⭐ Default" if sla.get("is_default") else ""
            active_badge = "✅ Active" if sla.get("is_active") else "❌ Inactive"
            
            with st.expander(f"{sla.get('name', 'Unnamed')} - {entity_display} {default_badge}"):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.write(f"**Entity Type**: {entity_display}")
                    st.write(f"**Max Hours**: {sla.get('max_review_hours', 'Unknown')}")
                    st.write(f"**Warning Threshold**: {sla.get('warning_threshold_pct', 'Unknown')}%")
                
                with col2:
                    st.write(f"**Status**: {active_badge}")
                    st.write(f"**Default**: {'Yes' if sla.get('is_default') else 'No'}")
                    escalation_count = len(sla.get("escalation_chain", []) or [])
                    st.write(f"**Escalation Levels**: {escalation_count}")
                
                with col3:
                    # Action buttons
                    col_edit, col_delete = st.columns(2)
                    
                    with col_edit:
                        if st.button("✏️ Edit", key=f"edit_sla_{sla['id']}"):
                            st.info("Edit functionality coming soon")
                    
                    with col_delete:
                        if st.button("🗑️ Delete", key=f"delete_sla_{sla['id']}"):
                            try:
                                with st.spinner("Deleting SLA rule..."):
                                    _api_delete(f"/api/v1/sla/rules/{sla['id']}")
                                st.success("✅ SLA rule deleted!")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Failed to delete SLA rule: {e}")
        
    except httpx.ConnectError:
        st.warning("⚠️ API unavailable. Showing demo data.")
        # Demo data
        st.write("**3 SLA Rules**")
        
        with st.expander("Global Default SLA - All Types ⭐ Default"):
            st.write("**Entity Type**: All Types")
            st.write("**Max Hours**: 120")
            st.write("**Warning Threshold**: 75%")
            st.write("**Status**: ✅ Active")
    except Exception as e:
        st.error(f"Failed to load SLA rules: {e}")


def render_sla_dashboard_page():
    """Render the main SLA Dashboard page."""
    user = check_authentication(AUTH_DISABLED)
    
    st.header("📈 SLA Dashboard")
    st.caption("Monitor and manage Service Level Agreements for review processes")
    
    # 4-tab layout
    tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Reviews", "Cycles", "Settings"])
    
    with tab1:
        render_overview_tab()
    
    with tab2:
        render_reviews_tab()
    
    with tab3:
        render_cycles_tab(user)
    
    with tab4:
        render_settings_tab(user)
