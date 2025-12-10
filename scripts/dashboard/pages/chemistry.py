"""Chemistry page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st

# Import models - use direct module import to avoid caching issues
import amprenta_rag.database.models as db_models

# Access models directly from module
Compound = db_models.Compound
HTSCampaign = db_models.HTSCampaign
BiochemicalResult = db_models.BiochemicalResult

from scripts.dashboard.db_session import db_session
from amprenta_rag.chemistry.compound_linking import (
    find_compounds_reversing_signature,
    find_signatures_affected_by_compound,
    link_compound_to_signature,
)
from amprenta_rag.chemistry.database import get_chemistry_db_path
import sqlite3


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
    tab1, tab2, tab3, tab4 = st.tabs(
        ["Compounds", "HTS Campaigns", "Biochemical Results", "Signature Links"]
    )

    with tab1:
        st.subheader("Compounds")
        with db_session() as db:
            search_term = st.text_input("Search compounds by ID or SMILES", "", key="compound_search")

            query = db.query(Compound)
            if search_term:
                query = query.filter(
                    (Compound.compound_id.ilike(f"%{search_term}%"))
                    | (Compound.smiles.ilike(f"%{search_term}%"))
                    | (Compound.inchi_key.ilike(f"%{search_term}%"))
                )

            compounds = query.order_by(Compound.created_at.desc()).limit(100).all()

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

        # Helper: fetch compounds from SQLite (chemistry DB)
        def _load_compound_ids() -> list[str]:
            try:
                conn = sqlite3.connect(str(get_chemistry_db_path()))
                rows = conn.execute("SELECT compound_id FROM compounds ORDER BY created_at DESC LIMIT 500").fetchall()
                return [r[0] for r in rows]
            except Exception:
                return []
            finally:
                try:
                    conn.close()
                except Exception:
                    pass

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
        st.markdown("### Find Compounds Reversing Signature")
        search_sig_id = st.text_input("Signature ID", key="siglink_search_signature")
        min_corr = st.slider(
            "Minimum correlation (negative)",
            min_value=-1.0,
            max_value=0.0,
            value=-0.5,
            step=0.05,
            key="siglink_min_corr",
        )
        if st.button("Search Reversing Compounds", key="siglink_search_btn"):
            if search_sig_id:
                try:
                    rows = find_compounds_reversing_signature(search_sig_id, min_correlation=min_corr)
                    if rows:
                        df = pd.DataFrame(rows, columns=["compound_id", "correlation", "p_value"])
                        df["effect_type"] = "reverses"
                        st.dataframe(df, hide_index=True, width='stretch')
                    else:
                        st.info("No reversing compounds found for this signature.")
                except Exception as e:
                    st.error(f"Error searching links: {e}")
            else:
                st.warning("Enter a signature ID to search.")

        st.markdown("---")
        st.markdown("### View Signature Links for Compound")
        compound_view_sel = st.selectbox(
            "Compound", compound_options, key="siglink_view_compound"
        )
        if st.button("Show Links", key="siglink_view_btn"):
            try:
                links = find_signatures_affected_by_compound(compound_view_sel)
                if links:
                    df_links = pd.DataFrame(
                        [
                            {
                                "signature_id": link.signature_id,
                                "effect_type": link.effect_type,
                                "correlation": link.correlation,
                                "p_value": link.p_value,
                                "evidence_source": link.evidence_source,
                            }
                            for link in links
                        ]
                    )
                    st.dataframe(df_links, hide_index=True, width='stretch')
                else:
                    st.info("No signature links for this compound.")
            except Exception as e:
                st.error(f"Error fetching links: {e}")
