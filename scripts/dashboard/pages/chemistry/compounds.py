"""Compounds, registration, and structure search tabs."""

from __future__ import annotations

import os

import pandas as pd
import httpx
import streamlit as st

from scripts.dashboard.components.mol3d_viewer import render_conformers_3d, render_molecule_3d
from scripts.dashboard.components.comment_widget import render_comments_widget

API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, json_body: dict, *, timeout: int = 120) -> dict:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    out = r.json()
    return out if isinstance(out, dict) else {"raw": out}


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
    elif endpoint == "logp":
        if 1 <= mean <= 3:
            level = "green"
        elif mean < 0 or mean > 5:
            level = "red"
        else:
            level = "yellow"
    return downgrade(level)


def _render_structural_alert_badge(light: str) -> None:
    tl = (light or "").upper()
    bg = {"GREEN": "#2E8B57", "YELLOW": "#F0AD4E", "RED": "#C73E1D"}.get(tl, "#888888")
    st.markdown(
        f"""
<div style="display:inline-block;padding:6px 12px;border-radius:999px;background:{bg};color:white;
font-weight:700;letter-spacing:0.5px;font-size:12px;">
{tl or "UNKNOWN"}
</div>
""",
        unsafe_allow_html=True,
    )


def render_compounds_tab(
    tab_compounds,
    tab_register,
    tab_structure,
    *,
    db_session,
    compound_model,
    register_compound,
    export_compounds,
    substructure_search,
    similarity_search,
    pharmacophore_search,
    get_pharmacophore_features,
) -> None:
    """Render Compounds, Register, and Structure Search tabs."""
    _render_compounds(tab_compounds, db_session, compound_model, export_compounds)
    _render_register(tab_register, register_compound)
    _render_structure_search(
        tab_structure,
        db_session=db_session,
        substructure_search=substructure_search,
        similarity_search=similarity_search,
        pharmacophore_search=pharmacophore_search,
        get_pharmacophore_features=get_pharmacophore_features,
        compound_model=compound_model,
    )


def _render_compounds(tab, db_session, compound_model, export_compounds):
    with tab:
        st.subheader("Compounds")
        search_term = st.text_input("Search compounds by ID or SMILES", "", key="compound_search")

        with db_session() as db:
            query = db.query(compound_model)
            if search_term:
                query = query.filter(
                    (compound_model.compound_id.ilike(f"%{search_term}%"))
                    | (compound_model.smiles.ilike(f"%{search_term}%"))
                )

            compounds = query.limit(100).all()

            st.metric("Compounds Found", len(compounds))
            if len(compounds) == 100:
                st.info("Showing first 100 compounds. Use search to narrow down.")

            if compounds:
                compound_data = []
                for compound in compounds:
                    hts_count = len(compound.hts_results) if compound.hts_results else 0
                    bio_count = len(compound.biochemical_results) if compound.biochemical_results else 0
                    compound_data.append(
                        {
                            "Compound ID": compound.compound_id,
                            "SMILES": compound.smiles[:50] + "..." if len(compound.smiles) > 50 else compound.smiles,
                            "InChI Key": compound.inchi_key or "-",
                            "Molecular Weight": compound.molecular_weight or "-",
                            "HTS Results": hts_count,
                            "Biochemical Results": bio_count,
                            "Created": compound.created_at.strftime("%Y-%m-%d"),
                        }
                    )
                df_compounds = pd.DataFrame(compound_data)
                st.dataframe(df_compounds, width="stretch", hide_index=True)

                st.markdown("### Export")
                col1, col2 = st.columns([1, 3])
                with col1:
                    export_format = st.selectbox("Format", ["csv", "json", "excel"], key="compound_export_format")
                with col2:
                    compound_ids = [c.compound_id for c in compounds]
                    mime_types = {
                        "csv": "text/csv",
                        "json": "application/json",
                        "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    }
                    file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                    try:
                        export_data = export_compounds(compound_ids, export_format, db)
                        st.download_button(
                            label=f"📥 Download Compounds ({export_format.upper()})",
                            data=export_data,
                            file_name=f"compounds.{file_extensions[export_format]}",
                            mime=mime_types[export_format],
                            key="export_compounds",
                        )
                    except Exception as e:
                        st.error(f"Export failed: {e}")

                st.markdown("---")
                st.subheader("Compound Details")

                for compound in compounds[:10]:
                    with st.expander(f"**{compound.compound_id}**"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**ID:** `{compound.id}`")
                            st.write(f"**Compound ID:** {compound.compound_id}")
                            st.write(f"**SMILES:** `{compound.smiles}`")
                            if compound.inchi_key:
                                st.write(f"**InChI Key:** `{compound.inchi_key}`")

                            # 3D viewer (lazy, inside its own expander)
                            with st.expander("3D View", expanded=False):
                                style = st.selectbox(
                                    "Style",
                                    options=["stick", "sphere", "line", "cartoon"],
                                    index=0,
                                    key=f"cmpd3d_style_{compound.id}",
                                )
                                show_confs = st.checkbox(
                                    "Show Conformers",
                                    value=False,
                                    key=f"cmpd3d_showconfs_{compound.id}",
                                )
                                if show_confs:
                                    render_conformers_3d(compound.smiles, n_conformers=5, style=style)
                                else:
                                    render_molecule_3d(compound.smiles, style=style)

                                pdb_cache_key = f"cmpd3d_pdb_{compound.id}"
                                if st.button("Generate PDB", key=f"cmpd3d_genpdb_{compound.id}"):
                                    try:
                                        out = _api_post(
                                            "/api/viz3d/conformers",
                                            {"smiles": compound.smiles, "n_conformers": 1, "optimize": True},
                                            timeout=120,
                                        )
                                        pdb = None
                                        if isinstance(out, dict):
                                            pdbs = out.get("pdb_strings") or []
                                            if isinstance(pdbs, list) and pdbs:
                                                pdb = str(pdbs[0])
                                        if not pdb:
                                            raise ValueError("No PDB returned")
                                        st.session_state[pdb_cache_key] = pdb
                                    except Exception as e:  # noqa: BLE001
                                        st.session_state.pop(pdb_cache_key, None)
                                        st.error(f"Failed to generate PDB: {e}")

                                pdb_str = st.session_state.get(pdb_cache_key)
                                if isinstance(pdb_str, str) and pdb_str.strip():
                                    st.download_button(
                                        "Download PDB",
                                        data=pdb_str.encode("utf-8"),
                                        file_name=f"{compound.compound_id}.pdb",
                                        mime="chemical/x-pdb",
                                        key=f"cmpd3d_dlpdb_{compound.id}",
                                    )
                            if compound.canonical_smiles:
                                st.write(f"**Canonical SMILES:** `{compound.canonical_smiles}`")
                        with col2:
                            st.write(f"**Created:** {compound.created_at.strftime('%Y-%m-%d %H:%M')}")
                            if compound.molecular_weight:
                                st.write(f"**Molecular Weight:** {compound.molecular_weight}")
                            if compound.molecular_formula:
                                st.write(f"**Molecular Formula:** {compound.molecular_formula}")
                            if compound.logp is not None:
                                st.write(f"**LogP:** {compound.logp}")
                            if compound.hbd_count is not None:
                                st.write(f"**HBD Count:** {compound.hbd_count}")
                            if compound.hba_count is not None:
                                st.write(f"**HBA Count:** {compound.hba_count}")

                            if st.button("View in Network", key=f"view_net_{compound.id}", type="secondary"):
                                st.session_state["ctn_prefill_compound_ids"] = [str(compound.id)]
                                st.session_state["selected_page"] = "Compound-Target Network"
                                st.rerun()

                        if st.button("Predict ADMET", key=f"admet_{compound.id}", type="secondary"):
                            try:
                                resp = _api_post(
                                    "/api/admet/predict",
                                    {"smiles": [compound.smiles], "endpoints": ["herg", "logs", "logp"], "include_uncertainty": True},
                                    timeout=120,
                                )
                                res0 = (resp.get("results") or [{}])[0]
                                preds = res0.get("predictions") or {}
                                rows = []
                                for ep, p in preds.items():
                                    if not isinstance(p, dict) or "mean" not in p:
                                        continue
                                    mean = float(p.get("mean"))
                                    std = p.get("std")
                                    light = _traffic_light(ep, mean, float(std) if std is not None else None)
                                    ci_low = p.get("ci_low")
                                    ci_high = p.get("ci_high")
                                    ci_str = "-"
                                    try:
                                        if ci_low is not None and ci_high is not None:
                                            ci_str = f"[{float(ci_low):.3f}, {float(ci_high):.3f}]"
                                    except Exception:  # noqa: BLE001
                                        ci_str = "-"
                                    rows.append(
                                        {
                                            "Endpoint": ep,
                                            "Mean": mean,
                                            "Std": std,
                                            "95% CI": ci_str,
                                            "In Domain": p.get("in_domain"),
                                            "Light": light,
                                        }
                                    )
                                if rows:
                                    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
                                else:
                                    st.info("No ADMET predictions returned (models may not be registered).")
                                st.markdown("[Open full ADMET Predictor](/?page=ADMET%20Predictor)")
                            except Exception as e:  # noqa: BLE001
                                st.error(f"ADMET prediction failed: {e}")

                        # Structural alerts (PAINS/Brenk/Lilly)
                        alert_cache_key = f"structural_alerts_check_{compound.id}"
                        c_btn, c_out = st.columns([1, 3])
                        with c_btn:
                            if st.button("Check Alerts", key=f"alerts_{compound.id}", type="secondary"):
                                with st.spinner("Checking structural alerts..."):
                                    try:
                                        out = _api_post("/api/alerts/check", {"smiles": compound.smiles, "filters": None}, timeout=60)
                                        st.session_state[alert_cache_key] = out
                                    except Exception as e:  # noqa: BLE001
                                        st.session_state[alert_cache_key] = {"smiles": compound.smiles, "error": str(e), "alerts": []}

                        with c_out:
                            cached = st.session_state.get(alert_cache_key)
                            if isinstance(cached, dict) and (cached.get("smiles") == compound.smiles):
                                if cached.get("error"):
                                    st.error(f"Alert check failed: {cached['error']}")
                                else:
                                    _render_structural_alert_badge(str(cached.get("traffic_light") or ""))
                                    summary = cached.get("summary") or {}
                                    if isinstance(summary, dict):
                                        st.caption(
                                            f"PAINS: {int(summary.get('pains_count', 0) or 0)} · "
                                            f"Brenk: {int(summary.get('brenk_count', 0) or 0)} · "
                                            f"Lilly: {int(summary.get('lilly_count', 0) or 0)}"
                                        )

                                    alerts = cached.get("alerts") or []
                                    if isinstance(alerts, list) and alerts:
                                        with st.expander(f"Alert details ({len(alerts)})", expanded=False):
                                            df_alerts = pd.DataFrame(alerts)
                                            if not df_alerts.empty:
                                                df_alerts = df_alerts.rename(
                                                    columns={
                                                        "alert_type": "Type",
                                                        "pattern_name": "Pattern Name",
                                                        "description": "Description",
                                                        "severity": "Severity",
                                                        "matched_smarts": "SMARTS",
                                                    }
                                                )
                                                cols = ["Type", "Pattern Name", "Description", "Severity", "SMARTS"]
                                                st.dataframe(
                                                    df_alerts[[c for c in cols if c in df_alerts.columns]],
                                                    hide_index=True,
                                                    use_container_width=True,
                                                )
                                    st.markdown("[Open full Structural Alerts](/?page=Structural%20Alerts)")

                        st.markdown(
                            f"[View in Graph Explorer](/?page=Graph%20Explorer&entity_type=compound&entity_id={compound.id})"
                        )

                        # Comments section
                        st.markdown("---")
                        render_comments_widget("compound", compound.id)
            else:
                st.info("No compounds found.")


def _render_register(tab, register_compound):
    with tab:
        st.subheader("Register Compound")
        
        # Add "Draw Structure" expander
        with st.expander("✏️ Draw Structure", expanded=False):
            from scripts.dashboard.components.ketcher_sketcher import render_ketcher_editor
            
            st.caption("Draw structure and click 'Get SMILES' to export")
            render_ketcher_editor(key="register_sketcher", width=700, height=500)
            st.info("Copy SMILES from editor above and paste into SMILES field below")
        
        name = st.text_input("Compound Name")
        smiles = st.text_input("SMILES")
        salt_form = st.selectbox(
            "Salt Form",
            ["None", "HCl", "sodium", "potassium", "free base", "other"],
            index=0,
        )
        if st.button("Register", type="primary"):
            if not name or not smiles:
                st.error("Name and SMILES are required.")
            else:
                try:
                    corp_id = register_compound(
                        name=name,
                        smiles=smiles,
                        salt_form=None if salt_form == "None" else salt_form,
                    )
                    st.success(f"Registered compound ID: {corp_id}")
                except Exception as e:
                    st.error(f"Registration failed: {e}")


def _render_structure_search(
    tab,
    *,
    db_session,
    substructure_search,
    similarity_search,
    pharmacophore_search,
    get_pharmacophore_features,
    compound_model,
):
    with tab:
        st.subheader("Structure Search")
        mode = st.radio("Search Mode", ["Substructure", "Similarity"], horizontal=True)
        query = st.text_input("SMARTS Pattern" if mode == "Substructure" else "Query SMILES")
        threshold = st.slider("Similarity Threshold", min_value=0.5, max_value=1.0, value=0.7, step=0.01)

        if st.button("Search", type="primary"):
            try:
                if mode == "Substructure":
                    results = substructure_search(query)
                else:
                    results = similarity_search(query, threshold=threshold)

                if not results:
                    st.info("No matches found.")
                else:
                    st.dataframe(results, hide_index=True, use_container_width=True)
            except Exception as e:
                st.error(f"Search failed: {e}")

        st.divider()
        st.markdown("### Pharmacophore Search")
        st.markdown("Search compounds by pharmacophore features (HBD, HBA, aromatic rings).")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("**Hydrogen Bond Donors (HBD)**")
            min_hbd = st.number_input("Min HBD", min_value=0, value=0, key="min_hbd")
            max_hbd = st.number_input("Max HBD", min_value=0, value=10, key="max_hbd")

        with col2:
            st.markdown("**Hydrogen Bond Acceptors (HBA)**")
            min_hba = st.number_input("Min HBA", min_value=0, value=0, key="min_hba")
            max_hba = st.number_input("Max HBA", min_value=0, value=10, key="max_hba")

        with col3:
            st.markdown("**Aromatic Rings**")
            min_aromatic = st.number_input("Min Aromatic", min_value=0, value=0, key="min_aromatic")
            max_aromatic = st.number_input("Max Aromatic", min_value=0, value=5, key="max_aromatic")

        if st.button("Search by Pharmacophore", type="primary", key="pharmacophore_search"):
            try:
                with db_session() as db:
                    results = pharmacophore_search(
                        min_hbd=min_hbd if min_hbd > 0 else None,
                        max_hbd=max_hbd if max_hbd < 10 else None,
                        min_hba=min_hba if min_hba > 0 else None,
                        max_hba=max_hba if max_hba < 10 else None,
                        min_aromatic=min_aromatic if min_aromatic > 0 else None,
                        max_aromatic=max_aromatic if max_aromatic < 5 else None,
                        db=db,
                        limit=100,
                    )

                    if not results:
                        st.info("No compounds found matching the criteria.")
                    else:
                        results_data = []
                        for compound in results:
                            features = get_pharmacophore_features(compound.smiles) if compound.smiles else {}
                            results_data.append(
                                {
                                    "Compound ID": compound.compound_id,
                                    "SMILES": compound.smiles[:50] if compound.smiles else "",
                                    "HBD": features.get("hbd", compound.hbd_count or 0),
                                    "HBA": features.get("hba", compound.hba_count or 0),
                                    "Aromatic Rings": features.get("aromatic_rings", 0),
                                    "LogP": compound.logp,
                                    "MW": compound.molecular_weight,
                                }
                            )

                        st.dataframe(pd.DataFrame(results_data), hide_index=True, use_container_width=True)
                        st.info(f"Found {len(results)} compounds")
            except Exception as e:
                st.error(f"Pharmacophore search failed: {e}")
                st.exception(e)

