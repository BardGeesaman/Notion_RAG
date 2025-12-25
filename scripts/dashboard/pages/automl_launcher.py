"""AutoML Templates launcher page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_automl_launcher_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("AutoML Templates")
    st.caption("Launch parameterized ML notebooks (open in JupyterHub or run headless via Papermill).")

    try:
        templates: List[Dict[str, Any]] = list(_api_get("/api/automl/templates") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load templates: {e}")
        return

    if not templates:
        st.info("No AutoML templates found.")
        return

    for tpl in templates:
        tid = str(tpl.get("id") or "")
        title = str(tpl.get("title") or tid)
        with st.expander(title):
            st.write(f"**Template ID:** `{tid}`")
            st.write(f"**Notebook path:** `{tpl.get('path')}`")
            st.link_button("Open in JupyterHub", str(tpl.get("jupyter_url") or ""))

            st.subheader("Parameters")
            dataset_id = st.text_input("dataset_id", value="", key=f"{tid}_dataset_id")
            target_column = st.text_input("target_column", value="target", key=f"{tid}_target")
            feature_columns = st.text_input(
                "feature_columns (comma-separated; blank = all numeric)", value="", key=f"{tid}_features"
            )
            test_size = st.slider("test_size", 0.1, 0.5, 0.2, 0.05, key=f"{tid}_test_size")
            n_clusters = st.slider("n_clusters (clustering)", 2, 20, 5, 1, key=f"{tid}_clusters")

            params: Dict[str, Any] = {"dataset_id": dataset_id}
            if "classification" in tid or "regression" in tid:
                params["target_column"] = target_column
                if feature_columns.strip():
                    params["feature_columns"] = [c.strip() for c in feature_columns.split(",") if c.strip()]
                params["test_size"] = float(test_size)
            if "clustering" in tid:
                if feature_columns.strip():
                    params["feature_columns"] = [c.strip() for c in feature_columns.split(",") if c.strip()]
                params["n_clusters"] = int(n_clusters)

            if st.button("Launch Training (Papermill)", key=f"{tid}_run", type="primary"):
                try:
                    out = _api_post("/api/automl/launch", {"template_id": tid, "params": params, "run_mode": "papermill"})
                    run_id = out.get("run_id")
                    if run_id:
                        st.success(f"Run created: {run_id}")
                        st.json(_api_get(f"/api/automl/runs/{run_id}"))
                    else:
                        st.warning("No run_id returned.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Launch failed: {e}")

    st.divider()
    st.subheader("Recent runs")
    try:
        runs = list(_api_get("/api/automl/runs") or [])
    except Exception:
        runs = []
    if runs:
        st.json(runs[-10:])
    else:
        st.caption("(none)")


__all__ = ["render_automl_launcher_page"]


