"""Dose-Response Explorer dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Tuple

import httpx
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _parse_pairs_text(text: str) -> Tuple[List[float], List[float]]:
    conc: List[float] = []
    resp: List[float] = []
    for raw in (text or "").splitlines():
        ln = raw.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = [p for p in ln.replace("\t", ",").replace(" ", ",").split(",") if p.strip()]
        if len(parts) < 2:
            continue
        try:
            conc.append(float(parts[0]))
            resp.append(float(parts[1]))
        except Exception:
            continue
    return conc, resp


def _parse_triplets_text(text: str) -> pd.DataFrame:
    rows = []
    for raw in (text or "").splitlines():
        ln = raw.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = [p for p in ln.replace("\t", ",").replace(" ", ",").split(",") if p.strip()]
        if len(parts) < 3:
            continue
        rows.append({"compound_id": parts[0], "concentration": parts[1], "response": parts[2]})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["concentration"] = pd.to_numeric(df["concentration"], errors="coerce")
    df["response"] = pd.to_numeric(df["response"], errors="coerce")
    df = df.dropna(subset=["compound_id", "concentration", "response"])
    return df


def _fit_curve_plot(
    conc: List[float],
    resp: List[float],
    fit: Dict[str, Any],
    *,
    title: str,
) -> go.Figure:
    curve_x = fit.get("curve_x") or []
    curve_y = fit.get("curve_y") or []
    ci_lo = fit.get("ci_lower")
    ci_hi = fit.get("ci_upper")

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=conc, y=resp, mode="markers", name="Data", marker=dict(size=8, color="#1f77b4")))
    if curve_x and curve_y:
        fig.add_trace(go.Scatter(x=curve_x, y=curve_y, mode="lines", name="Fit", line=dict(width=2, color="#d62728")))
    if curve_x and isinstance(ci_hi, list) and isinstance(ci_lo, list) and len(ci_hi) == len(curve_x) and len(ci_lo) == len(curve_x):
        fig.add_trace(go.Scatter(x=curve_x, y=ci_hi, mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip"))
        fig.add_trace(
            go.Scatter(
                x=curve_x,
                y=ci_lo,
                fill="tonexty",
                mode="lines",
                line=dict(width=0),
                fillcolor="rgba(214,39,40,0.2)",
                name="95% CI",
            )
        )
    fig.update_xaxes(type="log", title="Concentration (nM)")
    fig.update_yaxes(title="Response (%)")
    fig.update_layout(hovermode="x unified", template="plotly_white", title=title)
    return fig


def render_dose_response_explorer_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Dose-Response Explorer")
    st.caption("Fit dose-response curves and analyze time series via /api/explorer endpoints.")

    tab1, tab2, tab3, tab4 = st.tabs(["Single Curve", "Compare Curves", "Time Series", "Batch Analysis"])

    with tab1:
        st.subheader("Single Curve")
        text = st.text_area("Manual input (concentration,response per line)", height=140, key="dr_single_text")
        up = st.file_uploader("Or upload CSV (columns: concentration, response)", type=["csv"], key="dr_single_csv")

        model = st.selectbox("Model", options=["3PL", "4PL", "Bayesian 4PL"], index=1, key="dr_single_model")
        model_val = "bayesian_4pl" if model.lower().startswith("bayesian") else model

        conc: List[float] = []
        resp: List[float] = []
        if up is not None:
            try:
                df = pd.read_csv(up)
                if not {"concentration", "response"}.issubset(set(df.columns)):
                    st.error("CSV must contain columns: concentration, response")
                else:
                    conc = pd.to_numeric(df["concentration"], errors="coerce").dropna().astype(float).tolist()
                    resp = pd.to_numeric(df["response"], errors="coerce").dropna().astype(float).tolist()
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to read CSV: {e}")
        else:
            conc, resp = _parse_pairs_text(text)

        if st.button("Fit Curve", type="primary", key="dr_single_fit"):
            if not conc or not resp:
                st.warning("Provide concentration/response data first.")
            else:
                try:
                    out = _api_post(
                        "/api/explorer/dose-response/fit",
                        {"concentrations": conc, "responses": resp, "model": model_val},
                        timeout=180,
                    )
                    st.session_state["dr_single_fit"] = out
                except Exception as e:  # noqa: BLE001
                    st.session_state.pop("dr_single_fit", None)
                    st.error(f"Fit failed: {e}")

        fit = st.session_state.get("dr_single_fit")
        if isinstance(fit, dict):
            st.plotly_chart(_fit_curve_plot(conc, resp, fit, title="Dose-Response Fit"), use_container_width=True)
            params = {
                "EC50": fit.get("ec50"),
                "Hill": fit.get("hill_slope"),
                "Top": fit.get("top"),
                "Bottom": fit.get("bottom"),
                "R²": fit.get("r_squared"),
            }
            st.dataframe(pd.DataFrame([params]), hide_index=True, use_container_width=True)
            warns = fit.get("warnings") or []
            if warns:
                st.warning("\n".join([str(w) for w in warns]))
            st.download_button(
                "Download JSON",
                data=pd.Series(fit).to_json().encode("utf-8"),
                file_name="dose_response_fit.json",
                mime="application/json",
            )

    with tab2:
        st.subheader("Compare Curves")
        ups = st.file_uploader("Upload multiple CSVs (columns: concentration, response)", type=["csv"], accept_multiple_files=True, key="dr_cmp_csvs")
        txt2 = st.text_area("Or manual input (compound_id, concentration, response per line)", height=140, key="dr_cmp_text")

        fits: List[Dict[str, Any]] = []
        labels: List[str] = []

        if ups:
            for i, f in enumerate(ups[:50]):
                try:
                    df = pd.read_csv(f)
                    if not {"concentration", "response"}.issubset(set(df.columns)):
                        continue
                    c = pd.to_numeric(df["concentration"], errors="coerce").dropna().astype(float).tolist()
                    r = pd.to_numeric(df["response"], errors="coerce").dropna().astype(float).tolist()
                    fits.append({"concentrations": c, "responses": r, "model": "4PL"})
                    labels.append(getattr(f, "name", None) or f"curve_{i}")
                except Exception:
                    continue
        else:
            dfm = _parse_triplets_text(txt2)
            if not dfm.empty:
                for cid, g in dfm.groupby("compound_id"):
                    fits.append({"concentrations": g["concentration"].astype(float).tolist(), "responses": g["response"].astype(float).tolist(), "model": "4PL"})
                    labels.append(str(cid))

        if st.button("Compare", type="primary", key="dr_cmp_btn"):
            if not fits:
                st.warning("Provide at least one curve (CSV or manual).")
            else:
                try:
                    out = _api_post("/api/explorer/dose-response/compare", {"fits": fits[:50]}, timeout=300)
                    st.session_state["dr_cmp_out"] = out
                    st.session_state["dr_cmp_labels"] = labels[:50]
                except Exception as e:  # noqa: BLE001
                    st.session_state.pop("dr_cmp_out", None)
                    st.error(f"Compare failed: {e}")

        out = st.session_state.get("dr_cmp_out")
        labs = st.session_state.get("dr_cmp_labels") or []
        if isinstance(out, dict) and isinstance(out.get("results"), list):
            results = out["results"]
            fig = go.Figure()
            colors = px.colors.qualitative.Dark24
            for i, r in enumerate(results):
                if not isinstance(r, dict):
                    continue
                name = labs[i] if i < len(labs) else f"fit_{i}"
                cx = r.get("curve_x") or []
                cy = r.get("curve_y") or []
                if cx and cy:
                    fig.add_trace(go.Scatter(x=cx, y=cy, mode="lines", name=name, line=dict(width=2, color=colors[i % len(colors)])))
            fig.update_xaxes(type="log", title="Concentration (nM)")
            fig.update_yaxes(title="Response (%)")
            fig.update_layout(hovermode="x unified", template="plotly_white", title="Curve Comparison")
            st.plotly_chart(fig, use_container_width=True)

            rows = []
            for i, r in enumerate(results):
                if not isinstance(r, dict):
                    continue
                rows.append({"compound": labs[i] if i < len(labs) else f"fit_{i}", "ec50": r.get("ec50"), "hill": r.get("hill_slope"), "r2": r.get("r_squared")})
            df = pd.DataFrame(rows)
            st.dataframe(df, hide_index=True, use_container_width=True)
            st.download_button(
                "Download comparison CSV",
                data=df.to_csv(index=False).encode("utf-8"),
                file_name="dose_response_comparison.csv",
                mime="text/csv",
            )

    with tab3:
        st.subheader("Time Series")
        txt3 = st.text_area("Manual input (timepoint,value per line)", height=140, key="ts_text")
        up3 = st.file_uploader("Or upload CSV (columns: timepoint, value)", type=["csv"], key="ts_csv")

        smooth = st.selectbox("Smoothing", options=["None", "Savitzky-Golay", "LOWESS"], index=0, key="ts_smooth")
        detect_cp = st.checkbox("Detect Changepoints", value=False, key="ts_cp")

        tps: List[float] = []
        vals: List[float] = []
        if up3 is not None:
            try:
                df = pd.read_csv(up3)
                if not {"timepoint", "value"}.issubset(set(df.columns)):
                    st.error("CSV must contain columns: timepoint, value")
                else:
                    tps = pd.to_numeric(df["timepoint"], errors="coerce").dropna().astype(float).tolist()
                    vals = pd.to_numeric(df["value"], errors="coerce").dropna().astype(float).tolist()
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to read CSV: {e}")
        else:
            for raw in (txt3 or "").splitlines():
                ln = raw.strip()
                if not ln or ln.startswith("#"):
                    continue
                parts = [p for p in ln.replace("\t", ",").replace(" ", ",").split(",") if p.strip()]
                if len(parts) < 2:
                    continue
                try:
                    tps.append(float(parts[0]))
                    vals.append(float(parts[1]))
                except Exception:
                    continue

        if st.button("Analyze", type="primary", key="ts_analyze"):
            if not tps or not vals:
                st.warning("Provide timepoint/value data first.")
            else:
                method = None
                if smooth == "Savitzky-Golay":
                    method = "savgol"
                elif smooth == "LOWESS":
                    method = "lowess"
                try:
                    out = _api_post(
                        "/api/explorer/timeseries/analyze",
                        {
                            "values": vals,
                            "timepoints": tps,
                            "smooth_method": method,
                            "detect_changepoints": bool(detect_cp),
                            "changepoint_threshold": 2.0,
                        },
                        timeout=120,
                    )
                    st.session_state["ts_out"] = out
                except Exception as e:  # noqa: BLE001
                    st.session_state.pop("ts_out", None)
                    st.error(f"Analyze failed: {e}")

        out = st.session_state.get("ts_out")
        if isinstance(out, dict):
            slope = out.get("slope")
            pval = out.get("pvalue")
            direction = out.get("direction")
            c1, c2, c3 = st.columns(3)
            c1.metric("Slope", f"{slope:.4g}" if isinstance(slope, (int, float)) else "-")
            c2.metric("p-value", f"{pval:.4g}" if isinstance(pval, (int, float)) else "-")
            c3.metric("Direction", str(direction or "-"))

            fig = go.Figure()
            fig.add_trace(go.Scatter(x=tps, y=vals, mode="markers+lines", name="Raw", marker=dict(size=7, color="#1f77b4")))
            sm = out.get("smoothed_values")
            if isinstance(sm, list) and len(sm) == len(vals):
                fig.add_trace(go.Scatter(x=tps, y=sm, mode="lines", name="Smoothed", line=dict(width=2, color="#d62728")))
            cps = out.get("changepoint_indices")
            if isinstance(cps, list):
                for idx in cps:
                    try:
                        i = int(idx)
                        if 0 <= i < len(tps):
                            fig.add_vline(x=tps[i], line_width=1, line_dash="dash", line_color="#555555")
                    except Exception:
                        continue
            fig.update_layout(template="plotly_white", title="Time Series", hovermode="x unified")
            fig.update_xaxes(title="Timepoint")
            fig.update_yaxes(title="Value")
            st.plotly_chart(fig, use_container_width=True)

    with tab4:
        st.subheader("Batch Analysis")
        up4 = st.file_uploader(
            "Upload batch CSV (columns: compound_id, concentration_nm, response_percent)",
            type=["csv"],
            key="dr_batch_csv",
        )
        if st.button("Fit All", type="primary", key="dr_batch_fit"):
            if up4 is None:
                st.warning("Upload a batch CSV first.")
            else:
                try:
                    df = pd.read_csv(up4)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Failed to read CSV: {e}")
                    df = pd.DataFrame()
                req_cols = {"compound_id", "concentration_nm", "response_percent"}
                if not req_cols.issubset(set(df.columns)):
                    st.error("CSV must contain: compound_id, concentration_nm, response_percent")
                else:
                    fits = []
                    labels = []
                    for cid, g in df.groupby("compound_id"):
                        c = pd.to_numeric(g["concentration_nm"], errors="coerce").dropna().astype(float).tolist()
                        r = pd.to_numeric(g["response_percent"], errors="coerce").dropna().astype(float).tolist()
                        if len(c) < 4 or len(r) < 4:
                            continue
                        fits.append({"concentrations": c, "responses": r, "model": "4PL"})
                        labels.append(str(cid))
                        if len(fits) >= 50:
                            break
                    try:
                        out = _api_post("/api/explorer/dose-response/compare", {"fits": fits}, timeout=300)
                        st.session_state["dr_batch_out"] = out
                        st.session_state["dr_batch_labels"] = labels
                    except Exception as e:  # noqa: BLE001
                        st.session_state.pop("dr_batch_out", None)
                        st.error(f"Batch fit failed: {e}")

        bout = st.session_state.get("dr_batch_out")
        blabs = st.session_state.get("dr_batch_labels") or []
        if isinstance(bout, dict) and isinstance(bout.get("results"), list):
            rows = []
            for i, r in enumerate(bout["results"]):
                if not isinstance(r, dict):
                    continue
                warns = r.get("warnings") or []
                status = "ok" if not warns else "warning"
                rows.append(
                    {
                        "compound_id": blabs[i] if i < len(blabs) else f"fit_{i}",
                        "ec50": r.get("ec50"),
                        "hill": r.get("hill_slope"),
                        "r2": r.get("r_squared"),
                        "status": status,
                    }
                )
            df = pd.DataFrame(rows)
            st.dataframe(df, hide_index=True, use_container_width=True)
            st.download_button(
                "Download results CSV",
                data=df.to_csv(index=False).encode("utf-8"),
                file_name="dose_response_batch_results.csv",
                mime="text/csv",
            )

            # Heatmap (best-effort)
            try:
                df_hm = df.set_index("compound_id")[["ec50", "hill", "r2"]]
                fig = px.imshow(df_hm, aspect="auto", color_continuous_scale="Viridis", title="Parameters Heatmap")
                st.plotly_chart(fig, use_container_width=True)
            except Exception:
                pass


__all__ = ["render_dose_response_explorer_page"]


