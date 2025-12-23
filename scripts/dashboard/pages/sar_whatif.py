"""SAR What-If designer dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import plotly.express as px
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _parse_smiles_input(text: str) -> List[str]:
    out: List[str] = []
    for ln in (text or "").splitlines():
        t = ln.strip()
        if not t or t.startswith("#"):
            continue
        out.append(t)
    return out


def _risk_for_row(row: Dict[str, Any]) -> str:
    # Simple heuristic coloring
    logp = row.get("logp")
    logs = row.get("logs")
    herg = row.get("herg")

    red = False
    yellow = False

    try:
        if logp is not None and float(logp) >= 4.0:
            red = True
        elif logp is not None and float(logp) >= 3.0:
            yellow = True
    except Exception:
        pass

    try:
        if logs is not None and float(logs) <= -6.0:
            red = True
        elif logs is not None and float(logs) <= -4.0:
            yellow = True
    except Exception:
        pass

    try:
        if herg is not None and float(herg) >= 0.5:
            red = True
        elif herg is not None and float(herg) >= 0.3:
            yellow = True
    except Exception:
        pass

    if red:
        return "RED"
    if yellow:
        return "YELLOW"
    return "GREEN"


def render_sar_whatif_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("SAR What-If")
    st.caption("Predict ADMET properties for hypothetical analogs and explore scaffold hops.")

    tab1, tab2 = st.tabs(["Property Prediction", "Scaffold Hopping"])

    with tab1:
        st.subheader("Property Prediction")
        smiles_text = st.text_area(
            "SMILES (one per line)",
            height=160,
            placeholder="CCO\nc1ccccc1\nCC(=O)Oc1ccccc1C(=O)O",
        )
        if st.button("Predict", type="primary"):
            smiles_list = _parse_smiles_input(smiles_text)
            if not smiles_list:
                st.warning("Provide at least one SMILES.")
            else:
                try:
                    rows = _api_post("/api/v1/sar/predict", {"smiles_list": smiles_list}, timeout=300)
                    st.session_state["sar_whatif_last_rows"] = rows
                except Exception as e:  # noqa: BLE001
                    st.error(f"Prediction failed: {e}")

        rows = st.session_state.get("sar_whatif_last_rows") or []
        if rows:
            df = pd.DataFrame(rows)
            if "structure_img" in df.columns:
                # Streamlit ImageColumn expects URL; data URI works.
                df["structure_img"] = df["structure_img"].apply(
                    lambda b64: f"data:image/png;base64,{b64}" if isinstance(b64, str) and b64 else None
                )
            df["risk"] = [ _risk_for_row(r) for r in rows ]

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "structure_img": st.column_config.ImageColumn("Structure", width="small"),
                    "risk": st.column_config.TextColumn("Risk", help="GREEN/YELLOW/RED heuristic"),
                },
            )

            st.subheader("Property comparison")
            numeric_cols = [c for c in ("logp", "logs", "herg") if c in df.columns]
            if numeric_cols:
                melted = df[["smiles"] + numeric_cols].melt(id_vars=["smiles"], var_name="property", value_name="value")
                fig = px.bar(melted, x="smiles", y="value", color="property", barmode="group", title="ADMET properties")
                st.plotly_chart(fig, use_container_width=True)

    with tab2:
        st.subheader("Scaffold Hopping")
        ref_smiles = st.text_input("Reference SMILES", placeholder="c1ccccc1")
        try:
            tx = list(_api_get("/api/v1/sar/transformations") or [])
        except Exception:
            tx = []

        options = [t.get("id") for t in tx if t.get("id")]
        labels = {t.get("id"): t.get("label") for t in tx if t.get("id")}
        choice = st.selectbox(
            "Transformation",
            options=options,
            format_func=lambda x: f"{x} - {labels.get(x) or ''}",
        )
        if choice and ref_smiles and st.button("Generate Variants", type="primary"):
            try:
                out = _api_post("/api/v1/sar/scaffold-hop", {"smiles": ref_smiles, "transformation": choice}, timeout=120)
                products = out.get("products") or []
                st.session_state["sar_whatif_variants"] = [ref_smiles] + list(products)
            except Exception as e:  # noqa: BLE001
                st.error(f"Scaffold hop failed: {e}")

        variants = st.session_state.get("sar_whatif_variants") or []
        if variants:
            st.caption(f"Variants: {len(variants)}")
            try:
                rows2 = _api_post("/api/v1/sar/predict", {"smiles_list": variants}, timeout=300)
                df2 = pd.DataFrame(rows2)
                if "structure_img" in df2.columns:
                    df2["structure_img"] = df2["structure_img"].apply(
                        lambda b64: f"data:image/png;base64,{b64}" if isinstance(b64, str) and b64 else None
                    )
                df2["risk"] = [ _risk_for_row(r) for r in rows2 ]
                st.dataframe(
                    df2,
                    use_container_width=True,
                    hide_index=True,
                    column_config={"structure_img": st.column_config.ImageColumn("Structure", width="small")},
                )
            except Exception as e:  # noqa: BLE001
                st.error(f"Prediction failed: {e}")


__all__ = ["render_sar_whatif_page"]


