"""
Repository Health monitoring page.
"""

from __future__ import annotations

import datetime
import os
from typing import List, Dict

import httpx
import pandas as pd
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
SUMMARY_ENDPOINT = f"{API_BASE}/api/v1/catalog/summary"

HEALTH_COLORS = {
    "healthy": "#2ecc71",
    "warning": "#f1c40f",
    "stale": "#e74c3c",
}


# ---------------------------------------------------------------------------
# Data fetch
# ---------------------------------------------------------------------------
@st.cache_data(ttl=60)
def fetch_summary() -> Dict:
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(SUMMARY_ENDPOINT)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load repository summary: {e}")
        return {"repositories": []}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def health_score(status: str) -> int:
    status = (status or "").lower()
    if status == "healthy":
        return 100
    if status == "warning":
        return 60
    return 30


def render_status_cards(repos: List[Dict]):
    if not repos:
        st.info("No repository summary available.")
        return
    cols = st.columns(max(1, len(repos)))
    for col, repo in zip(cols, repos):
        color = HEALTH_COLORS.get((repo.get("health_status") or "").lower(), "#999")
        with col:
            st.markdown(
                f"""
                <div style="border:1px solid #ddd; border-radius:8px; padding:12px;">
                  <div style="font-weight:700; font-size:16px;">{repo.get('name','Unknown')}</div>
                  <div style="color:{color}; font-weight:700; margin:4px 0;">{repo.get('health_status','unknown').upper()}</div>
                  <div>Datasets: {repo.get('dataset_count', 0)}</div>
                  <div>Last sync: {repo.get('last_sync_date') or '—'}</div>
                  <div>Health score: {health_score(repo.get('health_status'))}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )


def render_sync_history(repos: List[Dict]):
    st.subheader("Sync History (last 30 days)")
    # Placeholder: no API for sync events yet
    st.info("Sync history endpoint not available. Add endpoint to visualize events.")


def render_health_metrics(repos: List[Dict]):
    st.subheader("Health Metrics (derived)")
    if not repos:
        st.info("No data for metrics.")
        return

    # Compute simple freshness metric: average age of last_sync_date
    ages = []
    now = datetime.datetime.now(datetime.timezone.utc)
    for r in repos:
        dt = r.get("last_sync_date")
        if isinstance(dt, str):
            try:
                dt = datetime.datetime.fromisoformat(dt)
            except ValueError:
                dt = None
        if dt:
            ages.append((now - dt).total_seconds() / 86400.0)
    data = {
        "Avg time between syncs (days)": "n/a (no history)",
        "Success rate": "n/a (no history)",
        "Data freshness (avg days since last sync)": f"{(sum(ages)/len(ages)):.2f}" if ages else "n/a",
    }
    st.table(pd.DataFrame(list(data.items()), columns=["Metric", "Value"]))


def render_actions():
    st.subheader("Quick Actions")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Trigger manual sync (placeholder)"):
            st.info("Manual sync trigger not implemented yet.")
    with col2:
        if st.button("View sync logs (placeholder)"):
            st.info("Sync logs view not implemented yet.")


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Repository Health", layout="wide")
st.title("Repository Health")

with st.spinner("Loading repository status..."):
    summary = fetch_summary()

repos = summary.get("repositories", []) if isinstance(summary, dict) else []

render_status_cards(repos)
st.divider()
render_sync_history(repos)
st.divider()
render_health_metrics(repos)
st.divider()
render_actions()


