"""Variant interpretation dashboard page."""

from __future__ import annotations

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


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 600) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_variant_analysis_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Variant Analysis")
    st.caption("Ingest VEP TSV variant sets, browse variants, compute gene burden and pathway enrichment.")

    tab1, tab2, tab3, tab4 = st.tabs(["Variant Sets", "Variants", "Gene Burden", "Pathway Enrichment"])

    sets: List[Dict[str, Any]] = []
    try:
        sets = list(_api_get("/api/variants/sets") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load variant sets: {e}")

    set_options = [s.get("id") for s in sets if s.get("id")]
    selected = st.session_state.get("variants_selected_set")
    if selected not in set_options:
        selected = set_options[0] if set_options else None

    with tab1:
        st.subheader("Variant Sets")
        if sets:
            st.dataframe(pd.DataFrame(sets), use_container_width=True, hide_index=True)
        else:
            st.info("No variant sets found.")

        st.divider()
        st.subheader("Ingest VEP TSV")
        name = st.text_input("Set name", placeholder="My Variant Set")
        vep_path = st.text_input("VEP TSV path", placeholder="/absolute/path/to/vep_output.tsv")
        if st.button("Ingest", type="primary"):
            try:
                out = _api_post("/api/variants/sets", {"name": name, "vep_tsv_path": vep_path}, timeout=1800)
                st.success(f"Ingested: {out.get('id')}")
                st.session_state["variants_selected_set"] = out.get("id")
            except Exception as e:  # noqa: BLE001
                st.error(f"Ingest failed: {e}")

    with tab2:
        st.subheader("Variant Table")
        sid = st.selectbox("Variant set", options=set_options, index=set_options.index(selected) if selected in set_options else 0 if set_options else None)
        st.session_state["variants_selected_set"] = sid
        gene = st.text_input("Filter gene (exact)", placeholder="TP53")
        impact = st.text_input("Filter impact (exact)", placeholder="HIGH")
        significance = st.text_input("Filter clinical significance (contains)", placeholder="pathogenic")
        limit = st.slider("Limit", min_value=10, max_value=1000, value=100, step=10)
        if sid:
            try:
                q = f"/api/variants/sets/{sid}/variants?limit={int(limit)}"
                if gene:
                    q += f"&gene={gene}"
                if impact:
                    q += f"&impact={impact}"
                if significance:
                    q += f"&significance={significance}"
                rows = list(_api_get(q) or [])
                if not rows:
                    st.info("No variants found for current filters.")
                else:
                    df = pd.DataFrame(rows)
                    st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load variants: {e}")
        else:
            st.info("No variant sets available.")

    with tab3:
        st.subheader("Gene Burden")
        sid = st.selectbox(
            "Variant set (burden)",
            options=set_options,
            index=set_options.index(selected) if selected in set_options else 0 if set_options else None,
            key="variants_burden_set",
        )
        min_score = st.slider("Min burden score", min_value=0.0, max_value=100.0, value=5.0, step=1.0)
        limit = st.slider("Top genes", min_value=10, max_value=200, value=50, step=10)
        if sid:
            try:
                rows = list(_api_get(f"/api/variants/sets/{sid}/burdens?min_score={float(min_score)}&limit={int(limit)}") or [])
                if not rows:
                    st.info("No gene burdens found (ingest a set first).")
                else:
                    df = pd.DataFrame(rows)
                    # Stacked bar: pathogenic/VUS/benign
                    df_plot = df.copy()
                    df_plot["gene"] = df_plot["gene_symbol"]
                    fig = px.bar(
                        df_plot,
                        x="gene",
                        y=["n_pathogenic", "n_vus", "n_benign"],
                        title="Gene burden counts (pathogenic/VUS/benign)",
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load gene burdens: {e}")
        else:
            st.info("No variant sets available.")

    with tab4:
        st.subheader("Pathway Enrichment")
        sid = st.selectbox(
            "Variant set (enrichment)",
            options=set_options,
            index=set_options.index(selected) if selected in set_options else 0 if set_options else None,
            key="variants_enrich_set",
        )
        min_b = st.slider("Min burden score for genes", min_value=0.0, max_value=100.0, value=5.0, step=1.0)
        if st.button("Run pathway enrichment"):
            if not sid:
                st.error("No variant set selected.")
            else:
                try:
                    rows = _api_post(
                        f"/api/variants/sets/{sid}/pathway-enrichment",
                        {"min_burden_score": float(min_b)},
                        timeout=120,
                    )
                    st.session_state["variants_last_enrichment"] = rows
                except Exception as e:  # noqa: BLE001
                    st.error(f"Enrichment failed: {e}")

        enr = st.session_state.get("variants_last_enrichment") or []
        if enr:
            df = pd.DataFrame(enr)
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("Run enrichment to see results.")


__all__ = ["render_variant_analysis_page"]


