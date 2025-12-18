"""
Alerts bell component for the dashboard.
"""

from __future__ import annotations

import os
import httpx
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
ALERTS_ENDPOINT = f"{API_BASE}/api/v1/alerts"


@st.cache_data(ttl=30)
def fetch_alerts(unread_only: bool = True):
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
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(f"{ALERTS_ENDPOINT}/read-all")
            resp.raise_for_status()
            fetch_alerts.clear()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to mark all alerts read: {e}")
        return False


def render_alerts_bell():
    """Render alerts bell in sidebar or header."""
    alerts = fetch_alerts(unread_only=True)
    count = len(alerts)
    if count > 0:
        st.markdown(f"🔔 **{count}**")
    else:
        st.markdown("🔔")

    with st.expander("Alerts", expanded=False):
        if not alerts:
            st.info("No unread alerts")
        else:
            for alert in alerts:
                col1, col2 = st.columns([4, 1])
                with col1:
                    st.write(f"New dataset: {alert.get('dataset_title') or alert.get('dataset_id')}")
                with col2:
                    if st.button("✓", key=f"read_{alert['id']}"):
                        if mark_alert_read(alert["id"]):
                            st.rerun()

            if st.button("Mark all read"):
                if mark_all_read():
                    st.rerun()

