"""
External Data Catalog dashboard page.

Displays repository summaries and a paginated table of external datasets.
"""

from __future__ import annotations

import os
import math
import io
import pandas as pd
import streamlit as st
import httpx


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
SUMMARY_ENDPOINT = f"{API_BASE}/api/v1/catalog/summary"
DATASETS_ENDPOINT = f"{API_BASE}/api/v1/catalog/datasets"
PAGE_SIZE = 50

HEALTH_COLORS = {
    "healthy": "#2ecc71",
    "warning": "#f1c40f",
    "stale": "#e74c3c",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@st.cache_data(ttl=60)
def fetch_summary():
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(SUMMARY_ENDPOINT)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load catalog summary: {e}")
        return {"repositories": []}


@st.cache_data(ttl=60)
def fetch_datasets(source: str | None, limit: int, offset: int):
    params = {"limit": limit, "offset": offset}
    if source and source != "All":
        params["source"] = source
    try:
        with httpx.Client(timeout=15) as client:
            resp = client.get(DATASETS_ENDPOINT, params=params)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load catalog datasets: {e}")
        return {"datasets": [], "total": 0}


def health_badge(health: str | None) -> str:
    color = HEALTH_COLORS.get((health or "").lower(), "#999")
    label = (health or "unknown").capitalize()
    return f"<span style='background:{color};color:#fff;padding:2px 8px;border-radius:8px;font-size:12px;'>{label}</span>"


def render_repo_cards(repos: list[dict]):
    if not repos:
        st.info("No repository summary available.")
        return
    cols = st.columns(max(1, len(repos)))
    for col, repo in zip(cols, repos):
        with col:
            badge = health_badge(repo.get("health_status"))
            last_sync = repo.get("last_sync_date") or "–"
            col.metric(label=repo.get("name", "Unknown"), value=int(repo.get("dataset_count", 0)))
            st.markdown(badge, unsafe_allow_html=True)
            st.caption(f"Last sync: {last_sync}")


def render_filters(sources: list[str]) -> dict:
    st.sidebar.header("Filters")
    src = st.sidebar.selectbox("Source", ["All"] + sources, index=0)
    search = st.sidebar.text_input("Search accession/title")
    refresh = st.sidebar.button("Refresh data")
    return {"source": src, "search": search, "refresh": refresh}


def apply_search(df: pd.DataFrame, term: str) -> pd.DataFrame:
    if not term:
        return df
    term_lower = term.lower()
    mask = (
        df["Accession"].fillna("").str.lower().str.contains(term_lower)
        | df["Title"].fillna("").str.lower().str.contains(term_lower)
    )
    return df[mask]


def render_table(datasets: list[dict], total: int, page: int, search: str):
    if not datasets:
        st.warning("No datasets found.")
        return

    df = pd.DataFrame(
        [
            {
                "Accession": d.get("accession"),
                "Title": d.get("title"),
                "Source": d.get("source"),
                "Import Date": d.get("created_at"),
                "Feature Count": d.get("feature_count", 0),
            }
            for d in datasets
        ]
    )

    df = apply_search(df, search)

    st.dataframe(df, use_container_width=True, hide_index=True)

    # Export
    csv_buf = io.StringIO()
    df.to_csv(csv_buf, index=False)
    st.download_button("Export CSV", csv_buf.getvalue(), file_name="external_catalog.csv", mime="text/csv")

    total_pages = max(1, math.ceil(total / PAGE_SIZE))
    st.caption(f"Page {page + 1} of {total_pages} — Total datasets: {total}")


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

st.set_page_config(page_title="External Catalog", layout="wide")
st.title("External Data Catalog")

# Fetch data
with st.spinner("Loading repository summary..."):
    summary = fetch_summary()
repos = summary.get("repositories", []) if isinstance(summary, dict) else []
render_repo_cards(repos)

sources = sorted({r.get("name") for r in repos if r.get("name")}) or []
filters = render_filters(sources)

# Simple page state
page = st.session_state.get("catalog_page", 0)
if filters["refresh"]:
    fetch_summary.clear()
    fetch_datasets.clear()
    page = 0

page = st.number_input("Page", min_value=1, value=page + 1, step=1) - 1
st.session_state["catalog_page"] = page

with st.spinner("Loading datasets..."):
    data = fetch_datasets(filters["source"], PAGE_SIZE, page * PAGE_SIZE)
datasets = data.get("datasets", []) if isinstance(data, dict) else []
total = data.get("total", len(datasets)) if isinstance(data, dict) else len(datasets)

render_table(datasets, total, page, filters["search"])

