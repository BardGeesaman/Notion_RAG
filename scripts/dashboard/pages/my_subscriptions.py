"""
My Subscriptions dashboard page.
"""

from __future__ import annotations

import json
import os
import io
import pandas as pd
import streamlit as st
import httpx

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
SUBS_ENDPOINT = f"{API_BASE}/api/v1/subscriptions"
PAGE_SIZE = 50

REPO_OPTIONS = ["GEO", "PRIDE", "MetaboLights", "MW", "all"]


# ---------------------------------------------------------------------------
# Fetchers with error handling
# ---------------------------------------------------------------------------
@st.cache_data(ttl=30)
def fetch_subscriptions():
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(SUBS_ENDPOINT)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load subscriptions: {e}")
        return []


def create_subscription(payload: dict):
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.post(SUBS_ENDPOINT, json=payload)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to create subscription: {e}")
        return None


def delete_subscription(sub_id: str):
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.delete(f"{SUBS_ENDPOINT}/{sub_id}")
            if resp.status_code not in (200, 204):
                resp.raise_for_status()
            return True
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to delete subscription: {e}")
        return False


def check_subscription(sub_id: str):
    try:
        with httpx.Client(timeout=15) as client:
            resp = client.post(f"{SUBS_ENDPOINT}/{sub_id}/check")
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to check subscription: {e}")
        return None


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------
def render_header():
    col1, col2 = st.columns([3, 1])
    with col1:
        st.title("My Subscriptions")
    with col2:
        st.button("New Subscription", key="new_sub_btn", type="primary", use_container_width=True, on_click=lambda: st.session_state.update(show_form=True))


def render_form():
    with st.expander("Create Subscription", expanded=st.session_state.get("show_form", False)):
        name = st.text_input("Name")
        source = st.selectbox("Repository source", REPO_OPTIONS, index=0)
        keywords = st.text_input("Keywords (comma-separated)")
        notify_email = st.checkbox("Notify by email", value=False)
        notify_in_app = st.checkbox("Notify in app", value=True)

        if st.button("Save", type="primary"):
            payload = {
                "name": name,
                "repository_source": source,
                "query_params": {"keywords": [k.strip() for k in keywords.split(",") if k.strip()]} if keywords else {},
                "notify_email": notify_email,
                "notify_in_app": notify_in_app,
            }
            created = create_subscription(payload)
            if created:
                st.success("Subscription created")
                st.session_state.show_form = False
                fetch_subscriptions.clear()
                st.rerun()


def render_table(subs: list[dict], page: int):
    if not subs:
        st.info("No subscriptions yet.")
        return

    df = pd.DataFrame(
        [
            {
                "Name": s.get("name"),
                "Source": s.get("repository_source"),
                "Keywords": ", ".join((s.get("query_params") or {}).get("keywords", [])),
                "Last Checked": s.get("last_checked"),
                "Status": "Active" if s.get("is_active") else "Paused",
                "ID": s.get("id"),
            }
            for s in subs
        ]
    )

    start = page * PAGE_SIZE
    end = start + PAGE_SIZE
    page_df = df.iloc[start:end]

    # Selection with radio buttons for clear mapping to actions
    options = [f"{row['Name']} ({row['Source']})" for _, row in page_df.iterrows()]
    ids = [row["ID"] for _, row in page_df.iterrows()]
    selection = st.radio("Select a subscription", options=options, index=0 if options else None)
    selected_id = ids[options.index(selection)] if options else None

    st.dataframe(
        page_df[["Name", "Source", "Keywords", "Last Checked", "Status"]],
        use_container_width=True,
        hide_index=True,
    )

    if selected_id:
        col1, col2 = st.columns([1, 1])
        with col1:
            if st.button("Check Now", key=f"check_{selected_id}"):
                with st.spinner("Checking..."):
                    res = check_subscription(selected_id)
                    if res is not None:
                        st.success("Check complete")
        with col2:
            confirm = st.checkbox("Confirm delete", key=f"confirm_{selected_id}")
            if st.button("Delete", key=f"del_{selected_id}", disabled=not confirm):
                with st.spinner("Deleting..."):
                    if delete_subscription(selected_id):
                        st.success("Deleted")
                        fetch_subscriptions.clear()
                        st.rerun()

    total_pages = (len(df) + PAGE_SIZE - 1) // PAGE_SIZE
    st.caption(f"Page {page+1} of {total_pages}")
    if total_pages > 1:
        new_page = st.number_input("Page", min_value=1, max_value=total_pages, value=page + 1, step=1) - 1
        if new_page != page:
            st.session_state.sub_page = new_page


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------
st.set_page_config(page_title="My Subscriptions", layout="wide")
render_header()
render_form()

page = st.session_state.get("sub_page", 0)

with st.spinner("Loading subscriptions..."):
    subs = fetch_subscriptions()

render_table(subs, page)

