"""Overview page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy import func

from amprenta_rag.database.models import (
    Compound,
    Dataset,
    Email,
    Experiment,
    Feature,
    HTSCampaign,
    Literature,
    Program,
    RAGChunk,
    Signature,
)
from scripts.dashboard.db_session import db_session


def render_overview_page() -> None:
    """
    Render the Overview page with statistics and recent datasets.

    Displays:
    - Overall statistics (datasets, programs, experiments, features, signatures)
    - Datasets by omics type distribution
    - Recent datasets list
    """
    st.header("📊 Overview")

    with db_session() as db:
        # Statistics - Core entities
        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            dataset_count = db.query(func.count(Dataset.id)).scalar()
            st.metric("Datasets", dataset_count)

        with col2:
            program_count = db.query(func.count(Program.id)).scalar()
            st.metric("Programs", program_count)

        with col3:
            experiment_count = db.query(func.count(Experiment.id)).scalar()
            st.metric("Experiments", experiment_count)

        with col4:
            feature_count = db.query(func.count(Feature.id)).scalar()
            st.metric("Features", feature_count)

        with col5:
            signature_count = db.query(func.count(Signature.id)).scalar()
            st.metric("Signatures", signature_count)

        st.markdown("---")

        # Statistics - Content entities
        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            literature_count = db.query(func.count(Literature.id)).scalar()
            st.metric("Literature", literature_count)

        with col2:
            email_count = db.query(func.count(Email.id)).scalar()
            st.metric("Emails", email_count)

        with col3:
            chunk_count = db.query(func.count(RAGChunk.id)).scalar()
            st.metric("RAG Chunks", chunk_count)

        with col4:
            compound_count = db.query(func.count(Compound.id)).scalar()
            st.metric("Compounds", compound_count)

        with col5:
            campaign_count = db.query(func.count(HTSCampaign.id)).scalar()
            st.metric("HTS Campaigns", campaign_count)

        # Datasets by Omics Type
        st.subheader("📈 Datasets by Omics Type")
        omics_counts = db.query(Dataset.omics_type, func.count(Dataset.id)).group_by(Dataset.omics_type).all()

        # Dataset creation over time
        dataset_timeline = (
            db.query(func.date(Dataset.created_at).label("date"), func.count(Dataset.id).label("count"))
            .group_by(func.date(Dataset.created_at))
            .order_by(func.date(Dataset.created_at))
            .all()
        )

    if omics_counts:
        df_omics = pd.DataFrame(omics_counts, columns=["Omics Type", "Count"])

        col1, col2 = st.columns(2)
        with col1:
            st.bar_chart(df_omics.set_index("Omics Type"))
        with col2:
            st.dataframe(df_omics, width='stretch', hide_index=True)

            # Export button
            csv = df_omics.to_csv(index=False)
            st.download_button(
                label="📥 Download CSV",
                data=csv,
                file_name="omics_distribution.csv",
                mime="text/csv",
            )

        # Dataset creation timeline
        if dataset_timeline:
            st.subheader("📅 Dataset Creation Timeline")
            df_timeline = pd.DataFrame([(str(d.date), d.count) for d in dataset_timeline], columns=["Date", "Count"])
            df_timeline["Date"] = pd.to_datetime(df_timeline["Date"])
            df_timeline = df_timeline.set_index("Date")
            st.line_chart(df_timeline)
    else:
        st.info("No datasets yet. Ingest some data to see statistics!")
        st.markdown(
            """
        **To get started:**
        1. Go to **Data Ingestion** page to upload datasets
        2. Or use CLI: `python scripts/ingest_lipidomics.py --file data.csv --create-page`
        3. Refresh this page to see your data
        """
        )

    # Recent Datasets
    st.subheader("Recent Datasets")
    with db_session() as db:
        recent_datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).limit(10).all()

        # Extract data while session is open
        dataset_data = []
        for ds in recent_datasets:
            dataset_data.append(
                {
                    "Name": ds.name,
                    "Omics Type": ds.omics_type,
                    "Created": ds.created_at.strftime("%Y-%m-%d %H:%M"),
                }
            )

    if dataset_data:
        df_recent = pd.DataFrame(dataset_data)
        st.dataframe(df_recent, width='stretch')
    else:
        st.info("No datasets found.")
