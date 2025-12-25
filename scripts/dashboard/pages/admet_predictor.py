"""ADMET Predictor dashboard page (uncertainty-aware)."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _parse_smiles(text: str, max_n: int = 100) -> List[str]:
    lines = [ln.strip() for ln in (text or "").splitlines()]
    smiles = [ln for ln in lines if ln]
    return smiles[: int(max_n)]


def _traffic_light(endpoint: str, mean: float, std: float | None) -> str:
    endpoint = (endpoint or "").lower()
    std_v = float(std) if std is not None else 0.0

    def downgrade(level: str) -> str:
        if std_v > 0.3:
            return {"green": "yellow", "yellow": "red", "red": "red"}.get(level, level)
        return level

    level = "yellow"
    if endpoint == "herg":
        if mean < 0.3:
            level = "green"
        elif mean > 0.7:
            level = "red"
    elif endpoint == "logs":
        if mean > -4:
            level = "green"
        elif mean < -6:
            level = "red"
        else:
            level = "yellow"
    elif endpoint == "logp":
        if 1 <= mean <= 3:
            level = "green"
        elif mean < 0 or mean > 5:
            level = "red"
        else:
            level = "yellow"

    return downgrade(level)


def _ci_str(ci_low: Any, ci_high: Any) -> str:
    try:
        lo = float(ci_low)
        hi = float(ci_high)
        return f"[{lo:.3f}, {hi:.3f}]"
    except Exception:  # noqa: BLE001
        return "-"


def render_admet_predictor_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("ADMET Predictor")
    st.caption("Uncertainty-aware ADMET predictions (ensemble + calibration + applicability domain).")

    tab1, tab2, tab3 = st.tabs(["Predict", "Calibration", "Model Info"])

    with tab1:
        st.subheader("Predict")
        smiles_text = st.text_area("SMILES (one per line, max 100)", height=140)
        show_uncertainty = st.checkbox("Show uncertainty", value=True)

        # Endpoint selection defaults
        endpoints_all = ["herg", "logs", "logp"]
        endpoints = st.multiselect("Endpoints", options=endpoints_all, default=endpoints_all)

        if st.button("Predict", type="primary"):
            smiles = _parse_smiles(smiles_text, max_n=100)
            if not smiles:
                st.error("Please enter at least one SMILES.")
            else:
                try:
                    out = _api_post(
                        "/api/admet/predict",
                        {"smiles": smiles, "endpoints": endpoints, "include_uncertainty": bool(show_uncertainty)},
                        timeout=120,
                    )
                    st.session_state["admet_predictor_last"] = out
                except Exception as e:  # noqa: BLE001
                    st.error(f"Prediction failed: {e}")

        out = st.session_state.get("admet_predictor_last")
        if isinstance(out, dict) and isinstance(out.get("results"), list):
            rows: List[Dict[str, Any]] = []
            errors: List[str] = []
            for r in out["results"]:
                if not isinstance(r, dict):
                    continue
                smi = str(r.get("smiles") or "")
                err = r.get("error")
                preds = r.get("predictions") or {}
                if err:
                    errors.append(str(err))
                    rows.append({"SMILES": smi, "Endpoint": "-", "Prediction": "-", "95% CI": "-", "Light": "red", "In Domain": False, "Error": err})
                    continue
                for ep, p in preds.items():
                    if not isinstance(p, dict):
                        continue
                    mean = p.get("mean")
                    std = p.get("std")
                    ci = _ci_str(p.get("ci_low"), p.get("ci_high")) if show_uncertainty else "-"
                    try:
                        mean_f = float(mean)
                    except Exception:  # noqa: BLE001
                        mean_f = float("nan")
                    std_f = float(std) if std is not None else None
                    light = _traffic_light(ep, mean_f, std_f)
                    rows.append(
                        {
                            "SMILES": (smi[:40] + "...") if len(smi) > 43 else smi,
                            "Endpoint": ep,
                            "Prediction": mean_f,
                            "95% CI": ci,
                            "Light": light,
                            "In Domain": p.get("in_domain"),
                            "Similarity": p.get("similarity"),
                            "Calibrated": p.get("calibrated"),
                            "Error": None,
                        }
                    )

            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.session_state["admet_predictor_df"] = df

            if errors:
                # Surface error text outside the dataframe so E2E can reliably assert it.
                st.error("; ".join(sorted(set(errors)))[:500])
        else:
            st.session_state.setdefault("admet_predictor_df", pd.DataFrame())

        # Always show export button (exports last results if available, else empty CSV).
        df_export = st.session_state.get("admet_predictor_df")
        if not isinstance(df_export, pd.DataFrame):
            df_export = pd.DataFrame()
        st.download_button(
            "Download CSV",
            data=df_export.to_csv(index=False).encode("utf-8"),
            file_name="admet_predictions.csv",
            mime="text/csv",
        )

    with tab2:
        st.subheader("Calibration")
        out = st.session_state.get("admet_predictor_last") or {}
        model_info = out.get("model_info") if isinstance(out, dict) else None
        if isinstance(model_info, dict):
            models = model_info.get("models") or {}
            metrics_rows: List[Dict[str, Any]] = []
            reliability_payloads: List[Dict[str, Any]] = []
            for ep, m in models.items():
                if not isinstance(m, dict):
                    continue
                metrics = m.get("metrics") or {}
                metrics_rows.append(
                    {
                        "endpoint": ep,
                        "ece": metrics.get("ece"),
                        "brier": metrics.get("brier"),
                        "auc": metrics.get("auc"),
                        "rmse": metrics.get("rmse"),
                        "r2": metrics.get("r2"),
                        "coverage_90": metrics.get("coverage_90"),
                    }
                )
                rd = metrics.get("reliability_diagram") or metrics.get("reliability")
                if isinstance(rd, dict):
                    reliability_payloads.append({"endpoint": ep, "rd": rd})
            if metrics_rows:
                st.dataframe(pd.DataFrame(metrics_rows), use_container_width=True, hide_index=True)
            else:
                st.info("No calibration metrics available yet (models may not be trained/registered).")

            # Reliability diagram visualization (if present)
            if reliability_payloads:
                try:
                    import plotly.graph_objects as go  # type: ignore
                except Exception:
                    st.info("Plotly not available; cannot render reliability diagram.")
                else:
                    for item in reliability_payloads:
                        ep = item["endpoint"]
                        rd = item["rd"]
                        bin_confs = rd.get("bin_confs") or []
                        bin_accs = rd.get("bin_accs") or []
                        fig = go.Figure()
                        fig.add_trace(go.Scatter(x=bin_confs, y=bin_accs, mode="lines+markers", name=f"{ep}"))
                        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", name="Perfect", line=dict(dash="dash")))
                        fig.update_layout(
                            title=f"Reliability Diagram ({ep})",
                            xaxis_title="Predicted probability",
                            yaxis_title="Observed frequency",
                            height=350,
                            margin=dict(l=40, r=20, t=50, b=40),
                        )
                        st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No calibration data available (reliability diagram missing from model metrics).")
        else:
            st.info("Run a prediction to view calibration metrics from the active model artifacts.")

        st.caption("Note: reliability diagrams are expected to be precomputed during training and stored in model metrics.")

    with tab3:
        st.subheader("Model Info")
        try:
            models = _api_get("/api/ml/models", timeout=60)
        except Exception as e:  # noqa: BLE001
            models = []
            st.error(f"Failed to load registry models: {e}")

        st.markdown("**Registered ADMET models (from MLModelRegistry)**")
        if isinstance(models, list) and models:
            df = pd.DataFrame([m for m in models if isinstance(m, dict) and str(m.get("name", "")).startswith("admet_")])
            if not df.empty:
                st.dataframe(df, use_container_width=True, hide_index=True)
            else:
                st.info("No admet_* models registered yet.")
        else:
            st.info("No models found.")

        last = st.session_state.get("admet_predictor_last")
        if isinstance(last, dict) and last.get("model_info"):
            st.markdown("**Active model info (from last prediction)**")
            st.json(last.get("model_info"))

        st.caption("Training script: scripts/train_admet_ensemble.py")


__all__ = ["render_admet_predictor_page"]


