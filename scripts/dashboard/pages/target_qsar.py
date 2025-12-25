"""Target QSAR dashboard page."""

from __future__ import annotations

import os
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


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def _parse_smiles(text: str, max_n: int = 100) -> List[str]:
    lines = [ln.strip() for ln in (text or "").splitlines()]
    smiles = [ln for ln in lines if ln]
    return smiles[: int(max_n)]


def _traffic_light(prob: float) -> str:
    try:
        p = float(prob)
    except Exception:  # noqa: BLE001
        return "YELLOW"
    if p > 0.7:
        return "GREEN"
    if p < 0.3:
        return "RED"
    return "YELLOW"


def render_target_qsar_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Target QSAR")
    st.caption("Predict compound activity (binary IC50 threshold) against per-target QSAR ensembles.")

    tab1, tab2, tab3 = st.tabs(["Available Models", "Predict", "Model Info"])

    # Shared: load targets list (cached in session_state)
    if st.button("Refresh Models", key="qsar_refresh_models"):
        st.session_state.pop("qsar_targets", None)

    targets = st.session_state.get("qsar_targets")
    if not isinstance(targets, list):
        try:
            targets = _api_get("/api/qsar/targets")
        except Exception as e:  # noqa: BLE001
            targets = []
            st.error(f"Failed to load QSAR targets: {e}")
        st.session_state["qsar_targets"] = targets

    target_names = [t.get("target") for t in targets if isinstance(t, dict) and t.get("target")]
    target_names = [str(x) for x in target_names]

    with tab1:
        st.subheader("Available Models")
        if not targets:
            st.info("No QSAR models found (or registry unavailable).")
        else:
            rows: List[Dict[str, Any]] = []
            for t in targets:
                if not isinstance(t, dict):
                    continue
                tgt = str(t.get("target") or "")
                metrics = t.get("metrics") if isinstance(t.get("metrics"), dict) else {}
                info = None
                try:
                    info_resp = _api_get(f"/api/qsar/targets/{tgt}/info")
                    info = info_resp.get("info") if isinstance(info_resp, dict) else None
                except Exception:
                    info = None
                rows.append(
                    {
                        "Target": tgt,
                        "Model Version": t.get("version"),
                        "AUC": (metrics or {}).get("auc"),
                        "Accuracy": (metrics or {}).get("accuracy"),
                        "Training Size": (info or {}).get("train_size") if isinstance(info, dict) else None,
                    }
                )
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)

    with tab2:
        st.subheader("Predict")
        if not target_names:
            st.info("No QSAR models available.")
        else:
            sel_targets = st.multiselect("Targets", options=target_names, default=target_names[:1])
            smi_text = st.text_area("SMILES (one per line, max 100)", height=160)

            if st.button("Predict", type="primary", key="qsar_predict_btn"):
                smiles = _parse_smiles(smi_text, max_n=100)
                if not smiles:
                    st.error("Please enter at least one SMILES.")
                elif not sel_targets:
                    st.error("Please select at least one target.")
                else:
                    try:
                        out = _api_post("/api/qsar/predict", {"smiles_list": smiles, "targets": sel_targets}, timeout=120)
                        st.session_state["qsar_predict_last"] = out
                    except Exception as e:  # noqa: BLE001
                        st.session_state["qsar_predict_last"] = None
                        st.error(f"QSAR prediction failed: {e}")

            out = st.session_state.get("qsar_predict_last")
            if isinstance(out, dict) and isinstance(out.get("results"), list):
                rows = []
                for r in out["results"]:
                    if not isinstance(r, dict):
                        continue
                    smi = str(r.get("smiles") or "")
                    err = r.get("error")
                    preds = r.get("predictions") or {}
                    if err:
                        rows.append({"SMILES": smi, "Target": "-", "Probability": None, "Active": None, "In-Domain": None, "Light": "RED", "Error": err})
                        continue
                    for tgt, p in preds.items():
                        if not isinstance(p, dict):
                            continue
                        prob = p.get("probability")
                        light = _traffic_light(prob)
                        rows.append(
                            {
                                "SMILES": (smi[:40] + "...") if len(smi) > 43 else smi,
                                "Target": tgt,
                                "Probability": prob,
                                "Active": p.get("active"),
                                "In-Domain": p.get("in_domain"),
                                "Light": light,
                                "Error": None,
                            }
                        )

                df = pd.DataFrame(rows)
                st.dataframe(df, use_container_width=True, hide_index=True)
                st.download_button(
                    "Download CSV",
                    data=df.to_csv(index=False).encode("utf-8"),
                    file_name="qsar_predictions.csv",
                    mime="text/csv",
                )

    with tab3:
        st.subheader("Model Info")
        if not target_names:
            st.info("No QSAR models available.")
        else:
            tgt = st.selectbox("Target", options=target_names, index=0)
            try:
                info = _api_get(f"/api/qsar/targets/{tgt}/info")
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load model info: {e}")
                info = None
            if isinstance(info, dict):
                st.json(info)


__all__ = ["render_target_qsar_page"]


