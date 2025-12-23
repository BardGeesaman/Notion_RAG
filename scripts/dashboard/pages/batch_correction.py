"""Batch correction dashboard page (ComBat)."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from scripts.dashboard.db_session import db_session


def render_batch_correction_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧪 Batch Effect Correction (ComBat)")
    st.caption("Select datasets, assign batch labels, and run ComBat correction.")

    with db_session() as db:
        datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).limit(200).all()

    if not datasets:
        st.info("No datasets available.")
        return

    ds_label = {f"{d.name} ({str(d.id)[:8]})": d for d in datasets}
    selected = st.multiselect("Select datasets", list(ds_label.keys()))
    if not selected:
        st.info("Select one or more datasets to continue.")
        return

    st.subheader("Assign batches")
    batch_map: Dict[str, Any] = {}
    for label in selected:
        ds = ds_label[label]
        batch = st.text_input(f"Batch label for {label}", value="batch1", key=f"batch_{ds.id}")
        batch_map[str(ds.id)] = batch

    method = st.selectbox("Method", ["combat"])

    api_url = os.environ.get("API_URL", "http://localhost:8000")
    endpoint = f"{api_url}/api/analysis/batch-correct"

    if st.button("Run batch correction", type="primary"):
        payload = {
            "datasets": [str(ds_label[label].id) for label in selected],
            "batch_map": batch_map,
            "method": method,
        }
        try:
            with httpx.Client(timeout=60) as client:
                resp = client.post(endpoint, json=payload)
            if resp.status_code >= 400:
                st.error(f"API error ({resp.status_code}): {resp.text}")
                return
            data = resp.json()
        except Exception as e:  # noqa: BLE001
            st.error(f"Request failed: {e}")
            return

        st.success("Batch correction completed.")
        st.json(data.get("stats") or {})

        preview = data.get("preview") or []
        if preview:
            st.subheader("Preview (first 20 rows)")
            st.dataframe(pd.DataFrame(preview), use_container_width=True, hide_index=True)

        corrected_id = data.get("corrected_dataset_id")
        if corrected_id:
            st.info(f"Created corrected dataset: {corrected_id}")


__all__ = ["render_batch_correction_page"]


