"""Single-Cell Viewer dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import httpx
import pandas as pd
import plotly.express as px
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_single_cell_viewer_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧫 Single-Cell Viewer")
    st.caption("Explore single-cell datasets: UMAP, clusters, and marker genes.")

    # Load datasets
    datasets: List[Dict[str, Any]] = []
    err: Optional[str] = None
    try:
        with httpx.Client(timeout=30) as client:
            r = client.get(f"{API_BASE}/api/single-cell/datasets", params={"limit": 200})
        r.raise_for_status()
        datasets = list(r.json() or [])
    except Exception as e:  # noqa: BLE001
        err = str(e)

    tab1, tab2, tab3, tab4 = st.tabs(["Dataset", "UMAP", "Clusters", "Marker Genes"])

    selected_id: Optional[str] = None
    with tab1:
        st.subheader("Dataset Selector")
        if err:
            st.error(f"Failed to load single-cell datasets: {err}")
        if not datasets:
            st.selectbox("Single-cell dataset", ["(none)"])
            st.info("Ingest a dataset via the API: POST /api/single-cell/ingest")
        else:
            opts = {}
            for d in datasets:
                sid = d.get("id")
                label = f"{sid[:8]} status={d.get('processing_status')} cells={d.get('n_cells')} genes={d.get('n_genes')}"
                opts[label] = sid
            label = st.selectbox("Single-cell dataset", list(opts.keys()))
            selected_id = opts[label]
            st.session_state["single_cell_selected_id"] = selected_id
            st.json(next((x for x in datasets if x.get("id") == selected_id), {}))

    selected_id = selected_id or st.session_state.get("single_cell_selected_id")

    with tab2:
        st.subheader("UMAP Viewer")
        if not selected_id:
            st.info("Select a dataset first.")
        else:
            gene_overlay = st.text_input("Gene expression overlay (optional)", placeholder="e.g., MS4A1")
            try:
                with httpx.Client(timeout=60) as client:
                    ur = client.get(f"{API_BASE}/api/single-cell/datasets/{selected_id}/umap", params={"limit": 20000})
                ur.raise_for_status()
                pts = list(ur.json() or [])
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load UMAP points: {e}")
                pts = []

            if not pts:
                st.info("No UMAP points available yet (dataset may still be processing).")
            else:
                df = pd.DataFrame(pts)
                color_col = "cluster_id"
                if gene_overlay.strip():
                    try:
                        with httpx.Client(timeout=60) as client:
                            er = client.get(
                                f"{API_BASE}/api/single-cell/datasets/{selected_id}/expression/{gene_overlay.strip()}"
                            )
                        if er.status_code == 200:
                            payload = er.json()
                            expr = pd.Series(payload.get("expression") or [])
                            obs_names = payload.get("obs_names") or []
                            # join by barcode/obs_names if lengths match
                            if len(expr) == len(obs_names):
                                expr_map = dict(zip(obs_names, expr))
                                df["expression"] = df["barcode"].map(expr_map)
                                color_col = "expression"
                    except Exception:
                        pass

                fig = px.scatter(
                    df,
                    x="umap_1",
                    y="umap_2",
                    color=color_col,
                    hover_data=["barcode", "cluster_id"],
                    height=650,
                )
                st.plotly_chart(fig, use_container_width=True)

    with tab3:
        st.subheader("Cluster Table")
        if not selected_id:
            st.info("Select a dataset first.")
        else:
            try:
                with httpx.Client(timeout=30) as client:
                    cr = client.get(f"{API_BASE}/api/single-cell/datasets/{selected_id}/clusters")
                cr.raise_for_status()
                rows = list(cr.json() or [])
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load clusters: {e}")
                rows = []
            if rows:
                st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            else:
                st.info("No clusters available yet.")

    with tab4:
        st.subheader("Marker Genes")
        if not selected_id:
            st.info("Select a dataset first.")
        else:
            cluster_filter = st.text_input("Filter by cluster_id (optional)", placeholder="e.g., 0")
            params: Dict[str, Any] = {"limit": 500}
            if cluster_filter.strip().isdigit():
                params["cluster_id"] = int(cluster_filter.strip())
            try:
                with httpx.Client(timeout=30) as client:
                    mr = client.get(f"{API_BASE}/api/single-cell/datasets/{selected_id}/markers", params=params)
                mr.raise_for_status()
                rows = list(mr.json() or [])
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load markers: {e}")
                rows = []
            if rows:
                df = pd.DataFrame(rows)
                cols = [c for c in ["gene_symbol", "cluster_id", "log2_fold_change", "pval_adj"] if c in df.columns]
                st.dataframe(df[cols], use_container_width=True, hide_index=True)
            else:
                st.info("No marker genes available yet.")


__all__ = ["render_single_cell_viewer_page"]


