"""Model Registry dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


def render_model_registry_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.set_page_config(page_title="Model Registry", layout="wide")
    st.title("🤖 ML Model Registry")

    api_url = os.environ.get("API_URL", "http://localhost:8000")
    base = f"{api_url}/api"

    models: List[Dict[str, Any]] = []
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(f"{base}/ml/models")
            if resp.ok:
                models = resp.json()
    except Exception:  # noqa: BLE001
        models = []

    if not models:
        st.info("No models registered yet. Train and register ADMET models to see them here.")
        return

    df = pd.DataFrame(models)

    col1, col2 = st.columns(2)
    with col1:
        type_opts = ["All"] + sorted(set(df["model_type"].dropna().tolist()))
        type_filter = st.selectbox("Model Type", type_opts)
    with col2:
        status_filter = st.selectbox("Status", ["All", "active", "archived", "training"])

    if type_filter != "All":
        df = df[df["model_type"] == type_filter]
    if status_filter != "All":
        df = df[df["status"] == status_filter]

    # Add monitoring links
    df_display = df.copy()
    df_display["monitoring"] = df_display["id"].apply(
        lambda model_id: f"[View Monitoring](?page=Model%20Monitoring&model_id={model_id})"
    )
    
    st.dataframe(
        df_display[["name", "version", "model_type", "framework", "status", "metrics", "monitoring"]],
        use_container_width=True,
        hide_index=True,
        column_config={
            "monitoring": st.column_config.LinkColumn("Monitoring", help="View model monitoring dashboard")
        }
    )

    if st.checkbox("Show model details"):
        selected_id = st.selectbox("Select model", df["id"].tolist())
        model_row = df[df["id"] == selected_id].iloc[0].to_dict()
        st.json(model_row)


__all__ = ["render_model_registry_page"]


