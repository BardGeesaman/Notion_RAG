"""ADMET Predictor dashboard page (uncertainty-aware)."""

from __future__ import annotations

import json
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


def render_shap_waterfall(shap_data: Dict[str, Any], title: str) -> None:
    """Render a SHAP waterfall chart from the Explain API payload."""
    try:
        import plotly.graph_objects as go  # type: ignore
    except Exception:
        st.info("Plotly not available; cannot render SHAP waterfall.")
        return

    top = shap_data.get("top_features") or []
    other_sum = float(shap_data.get("other_sum", 0.0))
    base_value = float(shap_data.get("base_value", 0.0))

    names = [str(f.get("name")) for f in top if isinstance(f, dict)]
    vals = [float(f.get("value", 0.0)) for f in top if isinstance(f, dict)]

    n_other = max(0, 2054 - len(names))
    names = names + [f"Other ({n_other} features)"]
    vals = vals + [other_sum]

    final_val = base_value + float(sum(vals))

    fig = go.Figure(
        go.Waterfall(
            x=["Base"] + names + ["Output"],
            measure=["absolute"] + ["relative"] * len(vals) + ["total"],
            y=[base_value] + vals + [final_val],
            connector={"line": {"color": "rgba(120,120,120,0.5)"}},
            increasing={"marker": {"color": "#2E8B57"}},  # green
            decreasing={"marker": {"color": "#C73E1D"}},  # red
            totals={"marker": {"color": "#444444"}},
        )
    )
    fig.update_layout(
        title=title,
        height=420,
        margin=dict(l=40, r=20, t=60, b=40),
        showlegend=False,
    )
    st.plotly_chart(fig, use_container_width=True)


def render_admet_predictor_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("ADMET Predictor")
    st.caption("Uncertainty-aware ADMET predictions (ensemble + calibration + applicability domain).")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(["Predict", "Calibration", "Model Info", "Explain", "Global Importance"])

    with tab1:
        st.subheader("Predict")
        smiles_text = st.text_area("SMILES (one per line, max 100)", height=140)
        show_uncertainty = st.checkbox("Show uncertainty", value=True)

        # Endpoint selection with categories
        st.write("**Select ADMET Endpoints:**")
        
        # Group endpoints by category for better UX
        toxicity_endpoints = ["herg"]
        metabolism_endpoints = ["cyp3a4", "cyp2d6", "cyp2c9"]
        solubility_endpoints = ["logs"]
        lipophilicity_endpoints = ["logp"]
        permeability_endpoints = ["bbb", "caco2"]
        clearance_endpoints = ["clearance"]
        
        endpoints_all = (
            toxicity_endpoints + metabolism_endpoints + solubility_endpoints + 
            lipophilicity_endpoints + permeability_endpoints + clearance_endpoints
        )
        
        # Create categorized display options
        endpoint_labels = {
            "herg": "hERG (Cardiotoxicity)",
            "cyp3a4": "CYP3A4 (Metabolism)", 
            "cyp2d6": "CYP2D6 (Metabolism)",
            "cyp2c9": "CYP2C9 (Metabolism)",
            "logs": "LogS (Solubility)",
            "logp": "LogP (Lipophilicity)", 
            "bbb": "BBB (Blood-Brain Barrier)",
            "caco2": "Caco-2 (Permeability)",
            "clearance": "Hepatic Clearance"
        }
        
        # Multi-select with descriptive labels
        endpoint_options = [f"{ep}: {endpoint_labels[ep]}" for ep in endpoints_all]
        selected_labels = st.multiselect(
            "Choose endpoints to predict",
            options=endpoint_options,
            default=endpoint_options[:3],  # Default to first 3 (herg, cyp3a4, cyp2d6)
            help="Select one or more ADMET properties to predict"
        )
        
        # Extract endpoint codes from selections
        endpoints = [label.split(":")[0] for label in selected_labels]

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

    with tab4:
        st.subheader("Explain")
        st.caption("Explain a single prediction using SHAP (top-K features + other).")

        smi = st.text_input("SMILES (single compound)", value="", max_chars=500)
        endpoint = st.selectbox("Endpoint", options=["herg", "logs", "logp"], index=0)
        top_k = st.slider("Top K features", min_value=5, max_value=20, value=10, step=1)

        if st.button("Explain", type="primary"):
            if not smi.strip():
                st.error("Please enter a SMILES string.")
            else:
                try:
                    out = _api_post(
                        "/api/admet/explain",
                        {"smiles": smi.strip(), "endpoint": endpoint, "top_k": int(top_k)},
                        timeout=120,
                    )
                    st.session_state["admet_explain_last"] = out
                except Exception as e:  # noqa: BLE001
                    # Keep UI usable even if API is unavailable.
                    st.session_state["admet_explain_last"] = {
                        "smiles": smi.strip(),
                        "endpoint": endpoint,
                        "prediction": {},
                        "shap": {"error": str(e), "top_features": [], "other_sum": 0.0, "base_value": 0.0},
                        "error": str(e),
                    }

        last = st.session_state.get("admet_explain_last") or {}
        if isinstance(last, dict) and last.get("smiles"):
            err = last.get("error")
            pred = last.get("prediction") if isinstance(last.get("prediction"), dict) else {}
            shap_payload = last.get("shap") if isinstance(last.get("shap"), dict) else None

            if err:
                st.error(str(err))

            if pred:
                c1, c2, c3 = st.columns(3)
                with c1:
                    try:
                        st.metric("Mean", f"{float(pred.get('mean')):.4f}")
                    except Exception:  # noqa: BLE001
                        st.metric("Mean", "-")
                with c2:
                    st.metric("95% CI", _ci_str(pred.get("ci_low"), pred.get("ci_high")))
                with c3:
                    st.metric("In Domain", str(pred.get("in_domain")))

            # Always render a stable SHAP section after an explain attempt, even if the API couldn't compute SHAP.
            if shap_payload is None:
                shap_payload = {"error": "No SHAP data available", "top_features": [], "other_sum": 0.0, "base_value": 0.0}

            if shap_payload.get("error"):
                st.warning(str(shap_payload.get("error")))

            render_shap_waterfall(shap_payload, title=f"SHAP Waterfall ({last.get('endpoint')})")

            st.download_button(
                "Download SHAP JSON",
                data=json.dumps(shap_payload, indent=2).encode("utf-8"),
                file_name="admet_shap.json",
                mime="application/json",
            )

            hist = st.session_state.setdefault("admet_shap_history", [])
            if isinstance(hist, list):
                hist.append({"endpoint": last.get("endpoint"), "top_features": shap_payload.get("top_features") or []})

    with tab5:
        st.subheader("Global Importance")
        st.caption("Requires multiple explanations to compute (MVP: accumulates top-K SHAP from Explain tab).")

        hist = st.session_state.get("admet_shap_history") or []
        if not isinstance(hist, list) or not hist:
            st.info("No SHAP data accumulated yet. Run a few explains first.")
        else:
            show_morgan = st.checkbox("Show Morgan bits", value=True)
            descriptors_only = st.checkbox("Show descriptors only", value=False)

            acc: Dict[str, float] = {}
            n = 0
            for item in hist:
                top = item.get("top_features") if isinstance(item, dict) else None
                if not isinstance(top, list):
                    continue
                n += 1
                for f in top:
                    if not isinstance(f, dict):
                        continue
                    name = str(f.get("name") or "")
                    val = float(f.get("value", 0.0))
                    acc[name] = acc.get(name, 0.0) + abs(val)

            if n == 0 or not acc:
                st.info("No usable SHAP top_features found yet.")
            else:
                imp = [{"name": k, "importance": v / float(n)} for k, v in acc.items()]

                if descriptors_only or not show_morgan:
                    imp = [x for x in imp if not str(x["name"]).startswith("MorganBit_")]

                imp.sort(key=lambda x: float(x["importance"]), reverse=True)
                top20 = imp[:20]

                if not top20:
                    st.info("No features to display under current filter.")
                else:
                    try:
                        import plotly.graph_objects as go  # type: ignore
                    except Exception:
                        st.info("Plotly not available; cannot render importance chart.")
                    else:
                        names = [x["name"] for x in top20][::-1]
                        vals = [float(x["importance"]) for x in top20][::-1]
                        fig = go.Figure(go.Bar(x=vals, y=names, orientation="h"))
                        fig.update_layout(
                            title="Top 20 features by mean |SHAP| (accumulated)",
                            height=520,
                            margin=dict(l=200, r=20, t=60, b=40),
                        )
                        st.plotly_chart(fig, use_container_width=True)


__all__ = ["render_admet_predictor_page"]


