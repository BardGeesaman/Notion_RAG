"""Sphingolipid Pathway Imbalance dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from scripts.dashboard.db_session import db_session


def render_sphingolipid_imbalance_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧬 Sphingolipid Pathway Imbalance")
    st.caption("Heuristic scoring for ceramide/sphingolipid pathway imbalance from dataset features.")

    with db_session() as db:
        datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).limit(200).all()

    if not datasets:
        st.info("No datasets available.")
        return

    ds_label = {f"{d.name} ({str(d.id)[:8]})": d for d in datasets}
    selected = st.multiselect("Select datasets", list(ds_label.keys()))
    pathway = st.selectbox("Pathway", ["ceramide"])

    if not selected:
        st.info("Select one or more datasets.")
        return

    api_url = os.environ.get("API_URL", "http://localhost:8000")
    endpoint = f"{api_url}/api/analysis/sphingolipid/imbalance"

    if st.button("Compute imbalance", type="primary"):
        payload = {"dataset_ids": [str(ds_label[label].id) for label in selected], "pathway": pathway}
        try:
            with httpx.Client(timeout=30) as client:
                resp = client.post(endpoint, json=payload)
            if resp.status_code >= 400:
                st.error(f"API error ({resp.status_code}): {resp.text}")
                return
            data: Dict[str, Any] = resp.json()
        except Exception as e:  # noqa: BLE001
            st.error(f"Request failed: {e}")
            return

        st.success(f"Score: {data.get('score'):.3f} ({data.get('method')})")
        st.json(data.get("stats") or {})

        matched = data.get("matched") or {}
        st.subheader("Matched signals")
        if isinstance(matched, dict):
            for k, v in matched.items():
                st.markdown(f"**{k}**")
                st.write(v)

        if "lipid_class_counts" in matched:
            st.subheader("Lipid class counts")
            st.dataframe(pd.DataFrame([matched["lipid_class_counts"]]), use_container_width=True, hide_index=True)


__all__ = ["render_sphingolipid_imbalance_page"]


