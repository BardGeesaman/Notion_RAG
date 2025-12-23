"""Multi-omics latent factor integration dashboard page."""

from __future__ import annotations

import json
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


def _safe_json(text: str) -> Any:
    if not text.strip():
        return None
    return json.loads(text)


def render_multi_omics_integration_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Multi-Omics Integration")
    st.caption("Create and run multi-omics latent factor experiments (MOFA-style).")

    tab1, tab2, tab3, tab4 = st.tabs(["Setup + Run", "Variance", "Factor Scatter", "Loadings"])

    experiments: List[Dict[str, Any]] = []
    try:
        experiments = list(_api_get("/api/multi-omics/experiments") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load experiments: {e}")
        experiments = []

    exp_options = [e.get("id") for e in experiments if e.get("id")]
    default_exp = st.session_state.get("multi_omics_selected_experiment")
    if default_exp not in exp_options:
        default_exp = exp_options[0] if exp_options else None

    with tab1:
        st.subheader("Existing experiments")
        if experiments:
            st.dataframe(pd.DataFrame(experiments), use_container_width=True, hide_index=True)
        else:
            st.info("No experiments found.")

        st.divider()
        st.subheader("Create experiment")
        name = st.text_input("Name", placeholder="My Multi-Omics Experiment")
        description = st.text_input("Description (optional)", placeholder="Short description")
        n_factors = st.number_input("n_factors", min_value=1, max_value=200, value=10, step=1)
        convergence_mode = st.selectbox("convergence_mode", options=["fast", "medium", "slow"], index=0)

        dataset_ids_text = st.text_area(
            "dataset_ids (JSON)",
            value='[{"omics_type":"transcriptomics","dataset_id":"<uuid>"},{"omics_type":"proteomics","dataset_id":"<uuid>"}]',
            height=120,
        )
        sample_mapping_text = st.text_area(
            "sample_mapping (JSON, optional)",
            value='[{"sample_id":"S1","transcriptomics":"RNA_A","proteomics":"PROT_A"}]',
            height=120,
        )

        c1, c2 = st.columns(2)
        with c1:
            if st.button("Create", type="primary"):
                try:
                    dataset_ids = _safe_json(dataset_ids_text)
                    sample_mapping = _safe_json(sample_mapping_text)
                    payload = {
                        "name": name,
                        "description": description or None,
                        "dataset_ids": dataset_ids,
                        "sample_mapping": sample_mapping,
                        "n_factors": int(n_factors),
                        "convergence_mode": str(convergence_mode),
                    }
                    out = _api_post("/api/multi-omics/experiments", payload, timeout=60)
                    st.success(f"Created experiment: {out.get('id')}")
                    st.session_state["multi_omics_selected_experiment"] = out.get("id")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Create failed: {e}")
        with c2:
            if st.button("Run selected"):
                eid = st.session_state.get("multi_omics_selected_experiment") or default_exp
                if not eid:
                    st.error("No experiment selected.")
                else:
                    try:
                        out = _api_post(f"/api/multi-omics/experiments/{eid}/run", {}, timeout=1800)
                        st.success(f"Run complete. Factors created: {out.get('factors_created')}")
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Run failed: {e}")

    with tab2:
        st.subheader("Variance explained per factor")
        eid = st.selectbox(
            "Experiment (variance)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_variance_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                if not factors:
                    st.info("No factors found (run experiment first).")
                else:
                    df = pd.DataFrame(factors)
                    df["total_variance"] = df["variance_explained"].apply(
                        lambda v: float(sum((v or {}).values())) if isinstance(v, dict) else None
                    )
                    fig = px.bar(df, x="factor_index", y="total_variance", title="Total variance explained (sum across views)")
                    st.plotly_chart(fig, use_container_width=True)
                    st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load variance: {e}")
        else:
            st.info("No experiments available.")

    with tab3:
        st.subheader("Factor scatter (Factor1 vs Factor2)")
        eid = st.selectbox(
            "Experiment (scatter)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_scatter_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                idxs = [int(f.get("factor_index")) for f in factors if f.get("factor_index") is not None]
                if len(idxs) < 2:
                    st.info("Need at least 2 factors (run experiment first).")
                else:
                    f1 = st.selectbox("X factor", options=idxs, index=0, key="multi_omics_f1")
                    f2 = st.selectbox("Y factor", options=idxs, index=1, key="multi_omics_f2")
                    s1 = pd.DataFrame(_api_get(f"/api/multi-omics/experiments/{eid}/factors/{f1}/scores") or [])
                    s2 = pd.DataFrame(_api_get(f"/api/multi-omics/experiments/{eid}/factors/{f2}/scores") or [])
                    if s1.empty or s2.empty:
                        st.info("No scores found (run experiment first).")
                    else:
                        s1 = s1.rename(columns={"score": "x"}).set_index("sample_id")
                        s2 = s2.rename(columns={"score": "y"}).set_index("sample_id")
                        merged = s1.join(s2, how="inner").reset_index()
                        fig = px.scatter(merged, x="x", y="y", hover_name="sample_id", title=f"Factor {f1} vs Factor {f2}")
                        st.plotly_chart(fig, use_container_width=True)
                        st.dataframe(merged, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load scatter: {e}")
        else:
            st.info("No experiments available.")

    with tab4:
        st.subheader("Top loadings")
        eid = st.selectbox(
            "Experiment (loadings)",
            options=exp_options,
            index=exp_options.index(default_exp) if default_exp in exp_options else 0 if exp_options else None,
            key="multi_omics_loadings_exp",
        )
        if eid:
            st.session_state["multi_omics_selected_experiment"] = eid
            try:
                factors = list(_api_get(f"/api/multi-omics/experiments/{eid}/factors") or [])
                idxs = [int(f.get("factor_index")) for f in factors if f.get("factor_index") is not None]
                if not idxs:
                    st.info("No factors found (run experiment first).")
                else:
                    fi = st.selectbox("Factor", options=idxs, index=0, key="multi_omics_load_factor")
                    limit = st.slider("Limit", min_value=10, max_value=200, value=50, step=10)
                    rows = list(
                        _api_get(f"/api/multi-omics/experiments/{eid}/factors/{fi}/loadings?limit={int(limit)}") or []
                    )
                    if not rows:
                        st.info("No loadings found.")
                    else:
                        df = pd.DataFrame(rows)
                        st.dataframe(df, use_container_width=True, hide_index=True)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load loadings: {e}")
        else:
            st.info("No experiments available.")


__all__ = ["render_multi_omics_integration_page"]


