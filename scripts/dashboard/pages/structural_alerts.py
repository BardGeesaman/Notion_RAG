"""Structural Alerts dashboard page (PAINS/Brenk/Lilly)."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def render_traffic_light(color: str) -> None:
    c = (color or "").upper()
    bg = {"GREEN": "#2E8B57", "YELLOW": "#F0AD4E", "RED": "#C73E1D"}.get(c, "#888888")
    st.markdown(
        f"""
<div style="display:inline-block;padding:6px 12px;border-radius:999px;background:{bg};color:white;
font-weight:700;letter-spacing:0.5px;font-size:12px;">
{c or "UNKNOWN"}
</div>
""",
        unsafe_allow_html=True,
    )


def _selected_filters(pains: bool, brenk: bool, lilly: bool) -> Optional[List[str]]:
    flt: List[str] = []
    if pains:
        flt.append("pains")
    if brenk:
        flt.append("brenk")
    if lilly:
        flt.append("lilly")
    return flt or None


def _alerts_df(alerts: List[Dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(alerts or [])
    if df.empty:
        return df
    ren = {
        "alert_type": "Type",
        "pattern_name": "Pattern Name",
        "description": "Description",
        "severity": "Severity",
        "matched_smarts": "SMARTS",
    }
    df = df.rename(columns=ren)
    cols = ["Type", "Pattern Name", "Description", "Severity", "SMARTS"]
    return df[[c for c in cols if c in df.columns]]


def _parse_batch_inputs(uploaded_file: Any, textarea: str) -> List[str]:
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to read CSV: {e}")
            return []
        cols = {str(c).strip().lower(): c for c in df.columns}
        if "smiles" not in cols:
            st.error("CSV must contain a 'smiles' column.")
            return []
        series = df[cols["smiles"]].dropna().astype(str)
        return [s.strip() for s in series.tolist() if s.strip()]

    lines = [ln.strip() for ln in (textarea or "").splitlines()]
    return [ln for ln in lines if ln]


def _download_csv_button(label: str, df: pd.DataFrame, filename: str) -> None:
    st.download_button(
        label,
        data=df.to_csv(index=False).encode("utf-8"),
        file_name=filename,
        mime="text/csv",
    )


def render_structural_alerts_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Structural Alerts")
    st.caption(
        "Run common medicinal-chemistry structural alert filters: "
        "PAINS (high-severity), Brenk (unwanted substructures), and a curated Lilly-style SMARTS list."
    )

    tab_single, tab_batch = st.tabs(["Single Compound", "Batch Check"])

    with tab_single:
        st.subheader("Single Compound")
        smi = st.text_input("SMILES", value="CCO", help="Enter a single SMILES string.")

        c1, c2, c3 = st.columns(3)
        with c1:
            pains = st.checkbox("PAINS", value=True)
        with c2:
            brenk = st.checkbox("Brenk", value=True)
        with c3:
            lilly = st.checkbox("Lilly", value=True)

        if st.button("Check Alerts", type="primary"):
            filters = _selected_filters(pains, brenk, lilly)
            try:
                out = _api_post("/api/alerts/check", {"smiles": smi, "filters": filters})
                st.session_state["structural_alerts_single"] = out
            except Exception as e:  # noqa: BLE001
                st.session_state["structural_alerts_single"] = None
                st.error(f"Alert check failed: {e}")

        out = st.session_state.get("structural_alerts_single")
        if isinstance(out, dict):
            render_traffic_light(str(out.get("traffic_light") or ""))

            summary = out.get("summary") or {}
            if isinstance(summary, dict):
                st.write(
                    {
                        "pains_count": int(summary.get("pains_count", 0) or 0),
                        "brenk_count": int(summary.get("brenk_count", 0) or 0),
                        "lilly_count": int(summary.get("lilly_count", 0) or 0),
                    }
                )

            if out.get("error"):
                st.error(str(out["error"]))

            alerts = out.get("alerts") or []
            with st.expander(f"Alert Details ({len(alerts)})", expanded=bool(alerts)):
                df = _alerts_df(alerts if isinstance(alerts, list) else [])
                if df.empty:
                    st.info("No alerts detected.")
                else:
                    st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_batch:
        st.subheader("Batch Check")
        st.caption("Upload a CSV with a `smiles` column, or paste one SMILES per line (max 100).")

        uploaded = st.file_uploader("Upload CSV", type=["csv"])
        text = st.text_area("SMILES list (one per line)", height=160)

        c1, c2, c3 = st.columns(3)
        with c1:
            pains = st.checkbox("PAINS (batch)", value=True)
        with c2:
            brenk = st.checkbox("Brenk (batch)", value=True)
        with c3:
            lilly = st.checkbox("Lilly (batch)", value=True)

        if st.button("Check Batch", type="primary"):
            smiles_list = _parse_batch_inputs(uploaded, text)
            smiles_list = [s for s in smiles_list if s][:100]
            if not smiles_list:
                st.error("Please provide at least one SMILES.")
            else:
                filters = _selected_filters(pains, brenk, lilly)
                try:
                    out = _api_post("/api/alerts/batch", {"smiles_list": smiles_list, "filters": filters}, timeout=120)
                    st.session_state["structural_alerts_batch"] = out
                except Exception as e:  # noqa: BLE001
                    st.session_state["structural_alerts_batch"] = None
                    st.error(f"Batch check failed: {e}")

        out = st.session_state.get("structural_alerts_batch")
        if isinstance(out, dict) and isinstance(out.get("results"), list):
            results = out["results"]
            df = pd.DataFrame(results)
            if not df.empty:
                # Flatten summary counts to columns for easier viewing.
                if "summary" in df.columns:
                    summ = df["summary"].apply(lambda x: x if isinstance(x, dict) else {})
                    df["pains_count"] = summ.apply(lambda d: int(d.get("pains_count", 0) or 0))
                    df["brenk_count"] = summ.apply(lambda d: int(d.get("brenk_count", 0) or 0))
                    df["lilly_count"] = summ.apply(lambda d: int(d.get("lilly_count", 0) or 0))
                st.dataframe(
                    df[["smiles", "traffic_light", "is_clean", "alert_count", "pains_count", "brenk_count", "lilly_count", "error"]]
                    if all(c in df.columns for c in ["smiles", "traffic_light", "is_clean", "alert_count"])
                    else df,
                    use_container_width=True,
                    hide_index=True,
                )

                clean_df = df[df.get("is_clean") == True][["smiles"]].copy()  # noqa: E712

                # One row per alert for flagged compounds.
                flagged_rows: List[Dict[str, Any]] = []
                for r in results:
                    if not isinstance(r, dict):
                        continue
                    if r.get("is_clean") is True:
                        continue
                    smi = r.get("smiles")
                    tl = r.get("traffic_light")
                    err = r.get("error")
                    alerts = r.get("alerts") or []
                    if not alerts:
                        flagged_rows.append(
                            {
                                "smiles": smi,
                                "traffic_light": tl,
                                "alert_type": "",
                                "pattern_name": "",
                                "description": "",
                                "severity": "",
                                "matched_smarts": "",
                                "error": err,
                            }
                        )
                    else:
                        for a in alerts:
                            if not isinstance(a, dict):
                                continue
                            flagged_rows.append(
                                {
                                    "smiles": smi,
                                    "traffic_light": tl,
                                    "alert_type": a.get("alert_type"),
                                    "pattern_name": a.get("pattern_name"),
                                    "description": a.get("description"),
                                    "severity": a.get("severity"),
                                    "matched_smarts": a.get("matched_smarts"),
                                    "error": err,
                                }
                            )
                flagged_df = pd.DataFrame(flagged_rows)

                c1, c2 = st.columns(2)
                with c1:
                    _download_csv_button("Download Clean Compounds", clean_df, "clean_compounds.csv")
                with c2:
                    if flagged_df.empty:
                        flagged_df = pd.DataFrame(
                            columns=[
                                "smiles",
                                "traffic_light",
                                "alert_type",
                                "pattern_name",
                                "description",
                                "severity",
                                "matched_smarts",
                                "error",
                            ]
                        )
                    _download_csv_button("Download Flagged Compounds", flagged_df, "flagged_compounds.csv")


__all__ = ["render_structural_alerts_page", "render_traffic_light"]


