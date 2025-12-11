"""Chemistry page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import plotly.express as px
import streamlit as st

# Import models - use direct module import to avoid caching issues
import amprenta_rag.database.models as db_models
from amprenta_rag.database.models import Compound as PgCompound
from amprenta_rag.database.base import get_db

# Access models directly from module
Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign
BiochemicalResult = db_models.BiochemicalResult
ADMEResult = db_models.ADMEResult
PKStudy = db_models.PKStudy
ToxicologyResult = db_models.ToxicologyResult

from scripts.dashboard.db_session import db_session
from amprenta_rag.chemistry.compound_linking import (
    link_compound_to_signature,
    get_compounds_for_signature,
)
from amprenta_rag.chemistry.registration import (
    check_duplicate,
    register_compound,
)
from amprenta_rag.chemistry.structure_search import substructure_search, similarity_search
from amprenta_rag.chemistry.sar_analysis import (
    get_compound_properties,
    get_activity_data,
    calculate_lipinski,
    detect_activity_cliffs,
)


def render_chemistry_page() -> None:
    """
    Render the Chemistry page with tabs for Compounds, HTS Campaigns, and Biochemical Results.

    Features:
    - Compounds: Search, filter, display properties
    - HTS Campaigns: Search, display results
    - Biochemical Results: Search, filter by assay
    - Export functionality (CSV)
    """
    st.header("⚗️ Chemistry")

    # Tabs for different chemistry entities
    tab1, tab_reg, tab_struct, tab_sar, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "Compounds",
            "Register Compound",
            "Structure Search",
            "SAR Analysis",
            "HTS Campaigns",
            "Biochemical Results",
            "Signature Links",
            "Lead Optimization",
        ]
    )

    with tab1:
        st.subheader("Compounds")
        search_term = st.text_input("Search compounds by ID or SMILES", "", key="compound_search")

        db_gen = get_db()
        db = next(db_gen)
        try:
            query = db.query(PgCompound)
            if search_term:
                query = query.filter(
                    (PgCompound.compound_id.ilike(f"%{search_term}%"))
                    | (PgCompound.smiles.ilike(f"%{search_term}%"))
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
                st.dataframe(df_compounds, width='stretch', hide_index=True)

                # Export button
                csv_compounds = df_compounds.to_csv(index=False)
                st.download_button(
                    label="📥 Download Compounds (CSV)",
                    data=csv_compounds,
                    file_name="compounds.csv",
                    mime="text/csv",
                    key="export_compounds",
                )

                st.markdown("---")
                st.subheader("Compound Details")

                for compound in compounds[:10]:  # Show first 10
                    with st.expander(f"**{compound.compound_id}**"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**ID:** `{compound.id}`")
                            st.write(f"**Compound ID:** {compound.compound_id}")
                            st.write(f"**SMILES:** `{compound.smiles}`")
                            if compound.inchi_key:
                                st.write(f"**InChI Key:** `{compound.inchi_key}`")
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
            else:
                st.info("No compounds found.")
        finally:
            db_gen.close()

    with tab_reg:
        st.subheader("Register Compound")
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

    with tab_struct:
        st.subheader("Structure Search")
        mode = st.radio("Search Mode", ["Substructure", "Similarity"], horizontal=True)
        if mode == "Substructure":
            query = st.text_input("SMARTS Pattern")
        else:
            query = st.text_input("Query SMILES")
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

    with tab_sar:
        st.subheader("SAR Analysis")
        props_df = get_compound_properties()
        if props_df.empty:
            st.info("No compounds available for SAR analysis.")
        else:
            # Lipinski summary
            lipinski_results = props_df.apply(
                lambda r: calculate_lipinski(r.get("smiles", "")),
                axis=1,
            )
            props_df["lipinski_compliant"] = [res.get("passes_ro5", False) for res in lipinski_results]
            compliant = int(props_df["lipinski_compliant"].sum())
            st.metric("Compounds", len(props_df))
            st.metric("Lipinski Compliant", compliant)

            col_a, col_b = st.columns(2)
            with col_a:
                fig_mw = px.histogram(props_df, x="molecular_weight", nbins=40, title="Molecular Weight")
                st.plotly_chart(fig_mw, use_container_width=True)
            with col_b:
                fig_logp = px.histogram(props_df, x="logp", nbins=40, title="LogP")
                st.plotly_chart(fig_logp, use_container_width=True)

            # Activity and cliffs
            act_df = get_activity_data()
            if not act_df.empty:
                merged = act_df.merge(
                    props_df[["compound_id", "smiles", "molecular_weight", "logp"]],
                    on="compound_id",
                    how="left",
                )
                merged_clean = merged.dropna(subset=["molecular_weight", "logp"])
                if not merged_clean.empty:
                    scatter = px.scatter(
                        merged_clean,
                        x="logp",
                        y="molecular_weight",
                        color="hit_flag",
                        title="LogP vs Molecular Weight (HTS)",
                        hover_data=["compound_id", "corporate_id"],
                    )
                    st.plotly_chart(scatter, use_container_width=True)

                cliff_thresh = st.slider(
                    "Cliff similarity threshold", min_value=0.5, max_value=1.0, value=0.8, step=0.01
                )
                cliffs = detect_activity_cliffs(merged, similarity_threshold=cliff_thresh)
                if cliffs:
                    st.markdown("**Activity Cliffs (top 50 by activity difference)**")
                    st.dataframe(cliffs[:50], hide_index=True, use_container_width=True)
                else:
                    st.info("No activity cliffs detected at the selected threshold.")
            else:
                st.info("No HTS activity data available for cliff detection.")

    with tab2:
        st.subheader("HTS Campaigns")
        with db_session() as db:
            search_term = st.text_input("Search campaigns by ID or name", "", key="campaign_search")

            query = db.query(HTSCampaign)
            if search_term:
                query = query.filter(
                    (HTSCampaign.campaign_id.ilike(f"%{search_term}%"))
                    | (HTSCampaign.campaign_name.ilike(f"%{search_term}%"))
                )

            campaigns = query.order_by(HTSCampaign.created_at.desc()).limit(100).all()

            st.metric("Campaigns Found", len(campaigns))
            if len(campaigns) == 100:
                st.info("Showing first 100 campaigns. Use search to narrow down.")

            if campaigns:
                campaign_data = []
                for campaign in campaigns:
                    result_count = len(campaign.results) if campaign.results else 0
                    campaign_data.append(
                        {
                            "Campaign ID": campaign.campaign_id,
                            "Name": (
                                campaign.campaign_name[:50] + "..."
                                if len(campaign.campaign_name) > 50
                                else campaign.campaign_name
                            ),
                            "Assay Type": campaign.assay_type or "-",
                            "Target": campaign.target or "-",
                            "Total Wells": campaign.total_wells or "-",
                            "Hit Count": campaign.hit_count or 0,
                            "Results": result_count,
                            "Created": campaign.created_at.strftime("%Y-%m-%d"),
                        }
                    )
                df_campaigns = pd.DataFrame(campaign_data)
                st.dataframe(df_campaigns, width='stretch', hide_index=True)

                # Export button
                csv_campaigns = df_campaigns.to_csv(index=False)
                st.download_button(
                    label="📥 Download Campaigns (CSV)",
                    data=csv_campaigns,
                    file_name="hts_campaigns.csv",
                    mime="text/csv",
                    key="export_campaigns",
                )

                st.markdown("---")
                st.subheader("Campaign Details")

                for campaign in campaigns[:10]:  # Show first 10
                    with st.expander(f"**{campaign.campaign_name}**"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**ID:** `{campaign.id}`")
                            st.write(f"**Campaign ID:** {campaign.campaign_id}")
                            st.write(f"**Name:** {campaign.campaign_name}")
                            if campaign.description:
                                st.write(f"**Description:** {campaign.description}")
                            if campaign.assay_type:
                                st.write(f"**Assay Type:** {campaign.assay_type}")
                            if campaign.target:
                                st.write(f"**Target:** {campaign.target}")
                        with col2:
                            st.write(f"**Created:** {campaign.created_at.strftime('%Y-%m-%d %H:%M')}")
                            if campaign.total_wells:
                                st.write(f"**Total Wells:** {campaign.total_wells}")
                            if campaign.hit_count:
                                st.write(f"**Hit Count:** {campaign.hit_count}")
                            if campaign.run_date:
                                st.write(f"**Run Date:** {campaign.run_date.strftime('%Y-%m-%d')}")

                            result_count = len(campaign.results) if campaign.results else 0
                            if result_count > 0:
                                st.write(f"**HTS Results:** {result_count}")
            else:
                st.info("No HTS campaigns found.")

    with tab3:
        st.subheader("Biochemical Results")
        with db_session() as db:
            search_term = st.text_input("Search by result ID or assay name", "", key="bio_search")

            query = db.query(BiochemicalResult)
            if search_term:
                query = query.filter(
                    (BiochemicalResult.result_id.ilike(f"%{search_term}%"))
                    | (BiochemicalResult.assay_name.ilike(f"%{search_term}%"))
                )

            results = query.order_by(BiochemicalResult.created_at.desc()).limit(100).all()

            st.metric("Results Found", len(results))
            if len(results) == 100:
                st.info("Showing first 100 results. Use search to narrow down.")

            if results:
                result_data = []
                for result in results:
                    result_data.append(
                        {
                            "Result ID": result.result_id,
                            "Assay Name": (
                                result.assay_name[:50] + "..." if len(result.assay_name) > 50 else result.assay_name
                            ),
                            "Target": result.target or "-",
                            "IC50": result.ic50 or "-",
                            "EC50": result.ec50 or "-",
                            "KI": result.ki or "-",
                            "KD": result.kd or "-",
                            "Created": result.created_at.strftime("%Y-%m-%d"),
                        }
                    )
                df_results = pd.DataFrame(result_data)
                st.dataframe(df_results, width='stretch', hide_index=True)

                # Export button
                csv_results = df_results.to_csv(index=False)
                st.download_button(
                    label="📥 Download Results (CSV)",
                    data=csv_results,
                    file_name="biochemical_results.csv",
                    mime="text/csv",
                    key="export_results",
                )

                st.markdown("---")
                st.subheader("Result Details")

                for result in results[:10]:  # Show first 10
                    with st.expander(f"**{result.assay_name}** - {result.result_id}"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**ID:** `{result.id}`")
                            st.write(f"**Result ID:** {result.result_id}")
                            st.write(f"**Assay Name:** {result.assay_name}")
                            if result.target:
                                st.write(f"**Target:** {result.target}")
                            if result.activity_type:
                                st.write(f"**Activity Type:** {result.activity_type}")
                        with col2:
                            st.write(f"**Created:** {result.created_at.strftime('%Y-%m-%d %H:%M')}")
                            if result.ic50:
                                st.write(f"**IC50:** {result.ic50} {result.units or ''}")
                            if result.ec50:
                                st.write(f"**EC50:** {result.ec50} {result.units or ''}")
                            if result.ki:
                                st.write(f"**KI:** {result.ki} {result.units or ''}")
                            if result.kd:
                                st.write(f"**KD:** {result.kd} {result.units or ''}")
                            if result.run_date:
                                st.write(f"**Run Date:** {result.run_date.strftime('%Y-%m-%d')}")
            else:
                st.info("No biochemical results found.")

    with tab4:
        st.subheader("Signature Links")

        # Helper: fetch compound IDs from PostgreSQL
        def _load_compound_ids() -> list[str]:
            db_gen = get_db()
            db = next(db_gen)
            try:
                compounds = db.query(PgCompound).order_by(PgCompound.created_at.desc()).limit(500).all()
                return [c.compound_id for c in compounds]
            except Exception:
                return []
            finally:
                db_gen.close()

        compound_options = _load_compound_ids()

        st.markdown("### Link Compound to Signature")
        col_l, col_r = st.columns(2)
        with col_l:
            compound_sel = st.selectbox("Compound", compound_options, key="siglink_compound_sel")
            signature_id = st.text_input("Signature ID (UUID)", key="siglink_signature_id")
            effect_type = st.selectbox(
                "Effect Type",
                ["reverses", "mimics", "partial", "unknown"],
                key="siglink_effect_type",
            )
        with col_r:
            correlation = st.number_input("Correlation", value=0.0, step=0.01, format="%.3f", key="siglink_corr")
            p_value = st.number_input("p-value", value=1.0, step=0.001, format="%.4f", key="siglink_pval")
            evidence_source = st.text_input("Evidence Source", key="siglink_evidence")

        if st.button("Create Link", key="siglink_create_btn"):
            if compound_sel and signature_id:
                try:
                    link_compound_to_signature(
                        compound_id=compound_sel,
                        signature_id=signature_id,
                        effect_type=effect_type,
                        correlation=correlation,
                        p_value=p_value,
                        evidence_source=evidence_source or None,
                    )
                    st.success(f"Linked compound {compound_sel} to signature {signature_id} ({effect_type}).")
                except Exception as e:
                    st.error(f"Error creating link: {e}")
            else:
                st.warning("Compound and Signature ID are required.")

        st.markdown("---")
        st.info("Use linking above; advanced signature searches removed.")

    with tab5:
        st.subheader("Lead Optimization")
        db_gen = get_db()
        db = next(db_gen)
        try:
            compound_choices = [
                c.compound_id for c in db.query(PgCompound).order_by(PgCompound.created_at.desc()).limit(500).all()
            ]
        except Exception:
            compound_choices = []
        finally:
            db_gen.close()

        col_lo1, col_lo2 = st.columns(2)

        with col_lo1:
            st.markdown("### ADME Data")
            adme_comp = st.selectbox("Compound (ADME)", compound_choices, key="lo_adme_comp")
            adme_assay = st.selectbox("Assay Type", ["permeability", "stability", "cyp_inhibition", "other"], key="lo_adme_assay")
            adme_value = st.number_input("Value", value=0.0, step=0.1, key="lo_adme_value")
            adme_unit = st.text_input("Unit", value="uM", key="lo_adme_unit")
            adme_conditions = st.text_area("Conditions (JSON)", key="lo_adme_conditions")
            if st.button("Save ADME", key="lo_adme_save"):
                db_gen = get_db()
                db = next(db_gen)
                try:
                    cond = None
                    if adme_conditions.strip():
                        import json
                        cond = json.loads(adme_conditions)
                    comp_obj = db.query(PgCompound).filter(PgCompound.compound_id == adme_comp).first() if adme_comp else None
                    rec = ADMEResult(
                        compound_id=comp_obj.id if comp_obj else None,
                        assay_type=adme_assay,
                        value=adme_value,
                        unit=adme_unit or None,
                        conditions=cond,
                    )
                    db.add(rec)
                    db.commit()
                    st.success("ADME entry saved.")
                except Exception as e:
                    db.rollback()
                    st.error(f"Failed to save ADME: {e}")
                finally:
                    db_gen.close()

        with col_lo2:
            st.markdown("### PK Studies")
            pk_comp = st.selectbox("Compound (PK)", compound_choices, key="lo_pk_comp")
            pk_species = st.text_input("Species", value="mouse", key="lo_pk_species")
            pk_route = st.text_input("Route", value="IV", key="lo_pk_route")
            pk_dose = st.number_input("Dose", value=0.0, step=0.1, key="lo_pk_dose")
            pk_dose_unit = st.text_input("Dose Unit", value="mg/kg", key="lo_pk_dose_unit")
            pk_auc = st.number_input("AUC", value=0.0, step=0.1, key="lo_pk_auc")
            pk_cmax = st.number_input("Cmax", value=0.0, step=0.1, key="lo_pk_cmax")
            pk_tmax = st.number_input("Tmax", value=0.0, step=0.1, key="lo_pk_tmax")
            pk_half = st.number_input("Half-life", value=0.0, step=0.1, key="lo_pk_half")
            pk_f = st.number_input("Bioavailability (%)", value=0.0, step=0.1, key="lo_pk_f")
            pk_cl = st.number_input("Clearance", value=0.0, step=0.1, key="lo_pk_cl")
            pk_vd = st.number_input("Volume of Distribution", value=0.0, step=0.1, key="lo_pk_vd")
            if st.button("Save PK", key="lo_pk_save"):
                db_gen = get_db()
                db = next(db_gen)
                try:
                    comp_obj = db.query(PgCompound).filter(PgCompound.compound_id == pk_comp).first() if pk_comp else None
                    rec = PKStudy(
                        compound_id=comp_obj.id if comp_obj else None,
                        species=pk_species or None,
                        route=pk_route or None,
                        dose=pk_dose or None,
                        dose_unit=pk_dose_unit or None,
                        auc=pk_auc or None,
                        c_max=pk_cmax or None,
                        t_max=pk_tmax or None,
                        half_life=pk_half or None,
                        bioavailability=pk_f or None,
                        clearance=pk_cl or None,
                        vd=pk_vd or None,
                    )
                    db.add(rec)
                    db.commit()
                    st.success("PK entry saved.")
                except Exception as e:
                    db.rollback()
                    st.error(f"Failed to save PK: {e}")
                finally:
                    db_gen.close()

        st.markdown("---")
        col_tox, col_profile = st.columns(2)

        with col_tox:
            st.markdown("### Toxicology")
            tox_comp = st.selectbox("Compound (Tox)", compound_choices, key="lo_tox_comp")
            tox_assay = st.text_input("Assay Type", value="herg", key="lo_tox_assay")
            tox_value = st.number_input("Value", value=0.0, step=0.1, key="lo_tox_value")
            tox_unit = st.text_input("Unit", value="uM", key="lo_tox_unit")
            tox_result = st.selectbox("Result", ["negative", "positive", "inconclusive"], key="lo_tox_result")
            tox_conditions = st.text_area("Conditions (JSON)", key="lo_tox_conditions")
            if st.button("Save Toxicology", key="lo_tox_save"):
                db_gen = get_db()
                db = next(db_gen)
                try:
                    cond = None
                    if tox_conditions.strip():
                        import json
                        cond = json.loads(tox_conditions)
                    comp_obj = db.query(PgCompound).filter(PgCompound.compound_id == tox_comp).first() if tox_comp else None
                    rec = ToxicologyResult(
                        compound_id=comp_obj.id if comp_obj else None,
                        assay_type=tox_assay or None,
                        value=tox_value or None,
                        unit=tox_unit or None,
                        result=tox_result or None,
                        conditions=cond,
                    )
                    db.add(rec)
                    db.commit()
                    st.success("Toxicology entry saved.")
                except Exception as e:
                    db.rollback()
                    st.error(f"Failed to save Toxicology: {e}")
                finally:
                    db_gen.close()

        with col_profile:
            st.markdown("### Compound Profile")
            profile_comp = st.selectbox("Compound (Profile)", compound_choices, key="lo_profile_comp")
            if profile_comp:
                db_gen = get_db()
                db = next(db_gen)
                try:
                    comp = db.query(PgCompound).filter(PgCompound.compound_id == profile_comp).first()
                    if comp:
                        adme = db.query(ADMEResult).filter(ADMEResult.compound_id == comp.id).all()
                        pk = db.query(PKStudy).filter(PKStudy.compound_id == comp.id).all()
                        tox = db.query(ToxicologyResult).filter(ToxicologyResult.compound_id == comp.id).all()

                        if adme:
                            st.markdown("**ADME Results**")
                            st.dataframe(
                                [
                                    {
                                        "assay_type": a.assay_type,
                                        "value": a.value,
                                        "unit": a.unit,
                                        "conditions": a.conditions,
                                    }
                                    for a in adme
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )
                        if pk:
                            st.markdown("**PK Studies**")
                            st.dataframe(
                                [
                                    {
                                        "species": p.species,
                                        "route": p.route,
                                        "auc": p.auc,
                                        "c_max": p.c_max,
                                        "t_max": p.t_max,
                                        "half_life": p.half_life,
                                        "bioavailability": p.bioavailability,
                                    }
                                    for p in pk
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )
                        if tox:
                            st.markdown("**Toxicology Results**")
                            st.dataframe(
                                [
                                    {
                                        "assay_type": t.assay_type,
                                        "value": t.value,
                                        "unit": t.unit,
                                        "result": t.result,
                                    }
                                    for t in tox
                                ],
                                hide_index=True,
                                use_container_width=True,
                            )

                        warnings_list = []
                        for t in tox:
                            if t.assay_type and t.assay_type.lower() == "herg" and t.value is not None and t.value < 10:
                                warnings_list.append("Potential hERG liability (IC50 < 10 uM)")
                            if t.assay_type and t.assay_type.lower() == "ames" and t.result and t.result.lower() == "positive":
                                warnings_list.append("Ames positive")
                        for p in pk:
                            if p.bioavailability is not None and p.bioavailability < 10:
                                warnings_list.append("Low bioavailability (<10%)")

                        if warnings_list:
                            for w in warnings_list:
                                st.error(w)
                        else:
                            st.success("No major liabilities detected.")
                finally:
                    db_gen.close()

