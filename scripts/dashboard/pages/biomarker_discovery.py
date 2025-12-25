"""Biomarker Discovery dashboard page."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def _try_read_matrix(file_path: str) -> pd.DataFrame:
    p = Path(str(file_path))
    if not p.exists():
        raise FileNotFoundError(f"Dataset file not found: {file_path}")

    for sep in (None, "\t", ",", ";"):
        try:
            df = pd.read_csv(p, sep=sep, engine="python")
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    raise ValueError(f"Unable to parse dataset file as a matrix: {file_path}")


def _detect_feature_column(df: pd.DataFrame) -> str:
    for cand in ("feature", "feature_id", "gene", "gene_id", "id", "name"):
        if cand in df.columns:
            return cand
    return str(df.columns[0])


def _extract_sample_ids_from_dataset_file(file_path: str) -> List[str]:
    raw = _try_read_matrix(file_path)
    feat_col = _detect_feature_column(raw)
    df2 = raw.copy().set_index(feat_col)
    num = df2.select_dtypes(include=["number"])
    if num.empty:
        num = df2.apply(pd.to_numeric, errors="coerce")
        num = num.dropna(axis=1, how="all")
    return [str(c) for c in num.columns.tolist()]


def render_biomarker_discovery_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Biomarker Discovery")
    st.caption("Run statistical + stability selection + CV importance and view a consensus feature ranking.")

    tab_setup, tab_methods, tab_results = st.tabs(["Setup", "Methods", "Results"])

    with tab_setup:
        st.subheader("Experiment + Sample Groups")

        if st.button("Refresh Experiments", key="bio_refresh_experiments"):
            st.session_state.pop("bio_experiments", None)

        exps = st.session_state.get("bio_experiments")
        if not isinstance(exps, list):
            try:
                exps = _api_get("/api/v1/experiments?limit=200")
            except Exception as e:  # noqa: BLE001
                exps = []
                st.error(f"Failed to load experiments: {e}")
            st.session_state["bio_experiments"] = exps

        exp_options: List[Dict[str, str]] = []
        for e in exps:
            if not isinstance(e, dict):
                continue
            eid = e.get("id")
            name = e.get("name") or "Unnamed"
            if eid:
                exp_options.append({"id": str(eid), "label": f"{name} ({str(eid)[:8]})"})

        if not exp_options:
            st.info("No experiments available.")
            return

        labels = [o["label"] for o in exp_options]
        sel_label = st.selectbox("Experiment", options=labels, index=0)
        sel_exp_id = next((o["id"] for o in exp_options if o["label"] == sel_label), None)
        st.session_state["bio_selected_experiment_id"] = sel_exp_id

        # Sample ID options: best-effort from first dataset file for the experiment.
        sample_options: List[str] = []
        if sel_exp_id:
            cache_key = f"bio_samples_{sel_exp_id}"
            cached = st.session_state.get(cache_key)
            if isinstance(cached, list):
                sample_options = [str(s) for s in cached]
            else:
                try:
                    dsets = _api_get(f"/api/v1/datasets?experiment_id={sel_exp_id}&limit=50")
                    file_path = None
                    if isinstance(dsets, list):
                        for ds in dsets:
                            if not isinstance(ds, dict):
                                continue
                            fps = ds.get("file_paths") or []
                            if isinstance(fps, list) and fps and fps[0]:
                                file_path = str(fps[0])
                                break
                    if file_path:
                        sample_options = _extract_sample_ids_from_dataset_file(file_path)
                    st.session_state[cache_key] = sample_options
                except Exception as e:  # noqa: BLE001
                    st.warning(f"Could not infer sample IDs from dataset files: {e}")
                    sample_options = []

        c1, c2 = st.columns(2)
        with c1:
            group1 = st.multiselect("Group 1 samples", options=sample_options, default=[])
        with c2:
            group2 = st.multiselect("Group 2 samples", options=sample_options, default=[])

        st.session_state["bio_group1"] = group1
        st.session_state["bio_group2"] = group2

        if not sample_options:
            st.info("If sample IDs could not be inferred, select them after fixing dataset file_paths, or paste IDs below.")
            g1_txt = st.text_area("Group 1 sample IDs (one per line)", height=90, key="bio_g1_txt")
            g2_txt = st.text_area("Group 2 sample IDs (one per line)", height=90, key="bio_g2_txt")
            if g1_txt.strip():
                st.session_state["bio_group1"] = [ln.strip() for ln in g1_txt.splitlines() if ln.strip()]
            if g2_txt.strip():
                st.session_state["bio_group2"] = [ln.strip() for ln in g2_txt.splitlines() if ln.strip()]

    with tab_methods:
        st.subheader("Methods + Run")
        try:
            methods = _api_get("/api/biomarker/methods")
        except Exception:
            methods = ["statistical", "stability", "importance"]

        default_methods = ["statistical", "stability", "importance"]
        sel_methods: List[str] = []
        for m in methods:
            m2 = str(m)
            checked = st.checkbox(m2, value=m2 in default_methods, key=f"bio_method_{m2}")
            if checked:
                sel_methods.append(m2)

        fdr = st.slider("FDR threshold", min_value=0.001, max_value=0.2, value=0.05, step=0.001)

        if st.button("Run Discovery", type="primary", key="bio_run_discovery"):
            exp_id = st.session_state.get("bio_selected_experiment_id")
            g1 = st.session_state.get("bio_group1") or []
            g2 = st.session_state.get("bio_group2") or []
            if not exp_id:
                st.error("Please select an experiment.")
            elif not g1 or not g2:
                st.error("Please select at least one sample in each group.")
            elif not sel_methods:
                st.error("Please select at least one method.")
            else:
                try:
                    out = _api_post(
                        "/api/biomarker/discover",
                        {
                            "experiment_id": exp_id,
                            "group1_samples": g1,
                            "group2_samples": g2,
                            "methods": sel_methods,
                            "fdr_threshold": float(fdr),
                        },
                        timeout=300,
                    )
                    st.session_state["bio_last_result"] = out
                except Exception as e:  # noqa: BLE001
                    st.session_state["bio_last_result"] = None
                    st.error(f"Biomarker discovery failed: {e}")

    with tab_results:
        st.subheader("Results")
        out = st.session_state.get("bio_last_result")
        if not isinstance(out, dict):
            st.info("Run discovery to see results.")
            return

        consensus = out.get("consensus_ranking") if isinstance(out.get("consensus_ranking"), list) else []
        if consensus:
            dfc = pd.DataFrame(consensus)
            st.markdown("**Consensus ranking** (average rank across selected methods)")
            st.dataframe(dfc, use_container_width=True, hide_index=True)
            st.download_button(
                "Download Consensus CSV",
                data=dfc.to_csv(index=False).encode("utf-8"),
                file_name="biomarker_consensus.csv",
                mime="text/csv",
            )
        else:
            st.info("No consensus ranking returned.")

        method_results = out.get("method_results") if isinstance(out.get("method_results"), dict) else {}
        if method_results:
            st.markdown("**Per-method results**")
            for name, rows in method_results.items():
                st.caption(str(name))
                if isinstance(rows, list):
                    dfm = pd.DataFrame(rows)
                    st.dataframe(dfm, use_container_width=True, hide_index=True)
                else:
                    st.json(rows)


__all__ = ["render_biomarker_discovery_page"]


