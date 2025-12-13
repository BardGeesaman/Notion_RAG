"""HTS Campaigns and Biochemical Results tabs."""

from __future__ import annotations

import pandas as pd
import streamlit as st


def render_assays_tab(
    tab_hts,
    tab_bio,
    *,
    db_session,
    hts_campaign_model,
    biochemical_result_model,
) -> None:
    """Render HTS Campaigns and Biochemical Results tabs."""
    _render_hts_campaigns(tab_hts, db_session, hts_campaign_model)
    _render_biochemical_results(tab_bio, db_session, biochemical_result_model)


def _render_hts_campaigns(tab, db_session, hts_campaign_model):
    with tab:
        st.subheader("HTS Campaigns")
        with db_session() as db:
            search_term = st.text_input("Search campaigns by ID or name", "", key="campaign_search")

            query = db.query(hts_campaign_model)
            if search_term:
                query = query.filter(
                    (hts_campaign_model.campaign_id.ilike(f"%{search_term}%"))
                    | (hts_campaign_model.campaign_name.ilike(f"%{search_term}%"))
                )

            campaigns = query.order_by(hts_campaign_model.created_at.desc()).limit(100).all()

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
                st.dataframe(df_campaigns, width="stretch", hide_index=True)

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

                for campaign in campaigns[:10]:
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


def _render_biochemical_results(tab, db_session, biochemical_result_model):
    with tab:
        st.subheader("Biochemical Results")
        with db_session() as db:
            search_term = st.text_input("Search by result ID or assay name", "", key="bio_search")

            query = db.query(biochemical_result_model)
            if search_term:
                query = query.filter(
                    (biochemical_result_model.result_id.ilike(f"%{search_term}%"))
                    | (biochemical_result_model.assay_name.ilike(f"%{search_term}%"))
                )

            results = query.order_by(biochemical_result_model.created_at.desc()).limit(100).all()

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
                st.dataframe(df_results, width="stretch", hide_index=True)

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

                for result in results[:10]:
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

