"""
Alerts bell component for the dashboard.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
import httpx
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
ALERTS_ENDPOINT = f"{API_BASE}/api/v1/alerts"
NOTIFICATIONS_ENDPOINT = f"{API_BASE}/api/v1/notifications"


@st.cache_data(ttl=30)
def fetch_alerts(unread_only: bool = True):
    """Fetch repository subscription alerts."""
    params = {"unread_only": str(unread_only).lower()}
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(ALERTS_ENDPOINT, params=params)
            resp.raise_for_status()
            alerts = resp.json()
            # Enrich for display
            for a in alerts:
                a.setdefault("dataset_title", a.get("dataset_id"))
            return alerts
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load alerts: {e}")
        return []


@st.cache_data(ttl=30)
def fetch_notifications(limit: int = 10):
    """Fetch activity notifications."""
    params = {"limit": limit}
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(NOTIFICATIONS_ENDPOINT, params=params)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load notifications: {e}")
        return []


@st.cache_data(ttl=30)
def fetch_notification_count():
    """Fetch unread notification count."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(f"{NOTIFICATIONS_ENDPOINT}/count")
            resp.raise_for_status()
            return resp.json().get("unread_count", 0)
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load notification count: {e}")
        return 0


def mark_alert_read(alert_id: str):
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(f"{ALERTS_ENDPOINT}/{alert_id}/read")
            resp.raise_for_status()
            fetch_alerts.clear()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to mark alert read: {e}")
        return False


def mark_all_read():
    """Mark all repository alerts as read."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(f"{ALERTS_ENDPOINT}/read-all")
            resp.raise_for_status()
            fetch_alerts.clear()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to mark all alerts read: {e}")
        return False


def mark_notification_read(notification_id: str):
    """Mark a single notification as read."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(f"{NOTIFICATIONS_ENDPOINT}/{notification_id}/read")
            resp.raise_for_status()
            fetch_notifications.clear()
            fetch_notification_count.clear()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to mark notification read: {e}")
        return False


def mark_all_notifications_read():
    """Mark all activity notifications as read."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(f"{NOTIFICATIONS_ENDPOINT}/read-all")
            resp.raise_for_status()
            fetch_notifications.clear()
            fetch_notification_count.clear()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to mark all notifications read: {e}")
        return False


def get_event_icon(event_type: str) -> str:
    """Get icon for activity event type."""
    icon_map = {
        "compound_added": "🧪",
        "experiment_created": "📊",
        "model_trained": "🤖",
        "hit_confirmed": "🎯",
        "status_changed": "📝",
        "notebook_reviewed": "📖",
    }
    return icon_map.get(event_type, "📋")


def format_time_ago(created_at: str) -> str:
    """Format time ago from ISO datetime string."""
    try:
        created_time = datetime.fromisoformat(created_at.replace('Z', '+00:00'))
        now = datetime.now(timezone.utc)
        diff = now - created_time
        
        if diff.days > 0:
            return f"{diff.days}d ago"
        elif diff.seconds > 3600:
            hours = diff.seconds // 3600
            return f"{hours}h ago"
        elif diff.seconds > 60:
            minutes = diff.seconds // 60
            return f"{minutes}m ago"
        else:
            return "Just now"
    except Exception:
        return "Unknown"


def render_alerts_bell():
    """Render alerts bell with activity and repository tabs."""
    # Get counts for badge
    alerts = fetch_alerts(unread_only=True)
    notification_count = fetch_notification_count()
    total_count = len(alerts) + notification_count
    
    # Bell icon with badge
    if total_count > 0:
        st.markdown(f"🔔 **{total_count}**")
    else:
        st.markdown("🔔")

    with st.expander("Notifications", expanded=False):
        tab1, tab2 = st.tabs(["Activity", "Repository"])
        
        with tab1:
            # Activity notifications tab
            notifications = fetch_notifications(limit=10)
            
            if not notifications:
                st.info("No recent activity")
            else:
                for notification in notifications:
                    activity_event = notification.get("activity_event", {})
                    event_type = activity_event.get("event_type", "unknown")
                    target_name = activity_event.get("target_name", "Unknown")
                    created_at = notification.get("created_at", "")
                    is_read = notification.get("is_read", True)
                    
                    # Create notification card
                    col1, col2, col3 = st.columns([0.5, 3, 1.5])
                    
                    with col1:
                        icon = get_event_icon(event_type)
                        st.markdown(f"**{icon}**")
                    
                    with col2:
                        # Title and description
                        title = f"{event_type.replace('_', ' ').title()}"
                        st.markdown(f"**{title}**")
                        st.caption(f"{target_name}")
                    
                    with col3:
                        time_ago = format_time_ago(created_at)
                        st.caption(time_ago)
                        if not is_read:
                            if st.button("✓", key=f"read_notif_{notification['id']}", help="Mark as read"):
                                if mark_notification_read(notification["id"]):
                                    st.rerun()
                
                # Mark all notifications read
                if notification_count > 0:
                    if st.button("Mark all activity read", use_container_width=True):
                        if mark_all_notifications_read():
                            st.rerun()
                
                # Link to full activity feed
                st.markdown("---")
                if st.button("📋 View Full Activity Feed", use_container_width=True):
                    st.switch_page("pages/activity_feed.py")
        
        with tab2:
            # Repository alerts tab (existing functionality)
            if not alerts:
                st.info("No unread repository alerts")
            else:
                for alert in alerts:
                    col1, col2, col3 = st.columns([4, 1, 1])
                    with col1:
                        ds_label = alert.get("dataset_title") or alert.get("dataset_id") or "Dataset"
                        if st.button(f"View {ds_label}", key=f"view_{alert['id']}", use_container_width=True):
                            if alert.get("dataset_id"):
                                # Auto-mark as read on navigation
                                mark_alert_read(alert["id"])
                                st.query_params["dataset_id"] = alert["dataset_id"]
                                st.switch_page("pages/dataset_details.py")
                    with col2:
                        if st.button("✓", key=f"read_alert_{alert['id']}"):
                            if mark_alert_read(alert["id"]):
                                st.rerun()
                    with col3:
                        st.caption("")

                # Mark all repository alerts read
                @st.dialog("Confirm")
                def confirm_mark_all_alerts():
                    st.write(f"Mark all {len(alerts)} repository alerts as read?")
                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button("Yes", use_container_width=True):
                            mark_all_read()
                            st.rerun()
                    with col2:
                        if st.button("Cancel", use_container_width=True):
                            st.rerun()

                if st.button("Mark all repository read"):
                    confirm_mark_all_alerts()

