"""CRISPR screen analysis dashboard page."""

from __future__ import annotations

import math
import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import plotly.express as px
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 300) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_crispr_analysis_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("CRISPR Analysis")
    st.caption("Create CRISPR screens, run MAGeCK, inspect hits and volcano plot.")

    tab1, tab2, tab3 = st.tabs(["Screens", "Results", "Volcano"])

    with tab1:
        st.subheader("Existing screens")
        screens: List[Dict[str, Any]] = []
        try:
            screens = list(_api_get("/api/crispr/screens") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load screens: {e}")

        if screens:
            st.dataframe(pd.DataFrame(screens), use_container_width=True, hide_index=True)
        else:
            st.info("No CRISPR screens found.")

        st.divider()
        st.subheader("Create screen")
        name = st.text_input("Name", placeholder="My CRISPR Screen")
        dataset_id = st.text_input("Dataset UUID", placeholder="UUID of dataset containing MAGeCK count matrix")
        library_type = st.text_input("Library type (optional)", placeholder="GeCKO v2 / Brunello / ...")
        cell_line = st.text_input("Cell line (optional)", placeholder="A375")
        treatment = st.text_input("Treatment (optional)", placeholder="Drug X")
        control_label = st.text_input("Control sample column name", placeholder="control")
        treatment_label = st.text_input("Treatment sample column name", placeholder="treatment")

        c1, c2 = st.columns(2)
        with c1:
            if st.button("Create", type="primary"):
                payload = {
                    "name": name,
                    "dataset_id": dataset_id,
                    "library_type": library_type or None,
                    "cell_line": cell_line or None,
                    "treatment": treatment or None,
                    "control_label": control_label,
                    "treatment_label": treatment_label,
                }
                try:
                    out = _api_post("/api/crispr/screens", payload, timeout=60)
                    st.success(f"Created screen: {out.get('id')}")
                    st.session_state["crispr_selected_screen"] = out.get("id")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Create failed: {e}")
        with c2:
            if st.button("Run analysis for selected"):
                sid = st.session_state.get("crispr_selected_screen")
                if not sid:
                    st.error("Select a screen first (Results tab selector).")
                else:
                    try:
                        out = _api_post(f"/api/crispr/screens/{sid}/analyze", {"method": "test"}, timeout=600)
                        st.success(f"Analysis complete. Results: {out.get('results')}")
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Analysis failed: {e}")

    # Shared selector
    screens_for_select: List[Dict[str, Any]] = []
    try:
        screens_for_select = list(_api_get("/api/crispr/screens") or [])
    except Exception:
        screens_for_select = []
    options = [s.get("id") for s in screens_for_select if s.get("id")]
    default = st.session_state.get("crispr_selected_screen")
    if default not in options:
        default = options[0] if options else None

    with tab2:
        st.subheader("Results")
        sid = st.selectbox("Screen", options=options, index=options.index(default) if default in options else 0 if options else None)
        st.session_state["crispr_selected_screen"] = sid
        if not sid:
            st.info("No screens available.")
        else:
            fdr_thr = st.slider("FDR threshold", min_value=0.0, max_value=1.0, value=0.05, step=0.01)
            limit = st.number_input("Limit", min_value=10, max_value=1000, value=100, step=10)
            if st.button("Refresh results"):
                pass
            try:
                rows = list(_api_get(f"/api/crispr/screens/{sid}/results?fdr_threshold={fdr_thr}&limit={int(limit)}") or [])
                if not rows:
                    st.info("No results found (run analysis first).")
                else:
                    df = pd.DataFrame(rows)
                    cols = ["rank", "gene_symbol", "beta_score", "fdr", "is_hit", "neg_lfc", "pos_lfc"]
                    st.dataframe(df[[c for c in cols if c in df.columns]], use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load results: {e}")

    with tab3:
        st.subheader("Volcano")
        sid2 = st.selectbox("Screen (volcano)", options=options, index=options.index(default) if default in options else 0 if options else None, key="crispr_volcano_screen")
        if not sid2:
            st.info("No screens available.")
        else:
            try:
                payload = _api_get(f"/api/crispr/screens/{sid2}/volcano-data")
                pts = payload.get("points") or []
                if not pts:
                    st.info("No volcano points found (run analysis first).")
                else:
                    df = pd.DataFrame(pts)
                    # Prefer neglog10_fdr computed by API; compute if missing
                    if "neglog10_fdr" not in df.columns and "fdr" in df.columns:
                        df["neglog10_fdr"] = df["fdr"].apply(lambda v: (-math.log10(v)) if v and v > 0 else None)
                    fig = px.scatter(
                        df,
                        x="neg_lfc",
                        y="neglog10_fdr",
                        color="is_hit",
                        hover_name="gene",
                        title="CRISPR Volcano (x=neg_lfc, y=-log10(FDR))",
                    )
                    st.plotly_chart(fig, use_container_width=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load volcano data: {e}")


__all__ = ["render_crispr_analysis_page"]


