"""Overview page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy import func

from scripts.dashboard.utils.accessibility import (
    render_skip_link,
    add_heading_structure,
    ensure_minimum_contrast
)

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
from amprenta_rag.utils.activity import (
    get_activity_stats,
    get_recent_compounds,
    get_recent_discoveries,
    get_recent_experiments,
)
from amprenta_rag.utils.widgets import (
    get_sample_count,
    get_discovery_count,
)
from scripts.dashboard.db_session import db_session
from scripts.dashboard.utils.cache import (
    fetch_programs, 
    fetch_datasets, 
    fetch_experiments, 
    fetch_compounds,
    clear_all_caches
)


def render_overview_page() -> None:
    """
    Render the Overview page with statistics and recent datasets.

    Displays:
    - Overall statistics (datasets, programs, experiments, features, signatures) - CACHED
    - Datasets by omics type distribution
    - Recent datasets list
    - Cache refresh button
    """
    # Add accessibility features
    render_skip_link("main-overview-content")
    ensure_minimum_contrast()
    
    # Add main content landmark and heading
    st.markdown(
        """
        <main id="main-overview-content" role="main" aria-label="Platform overview dashboard">
        </main>
        """,
        unsafe_allow_html=True
    )
    
    add_heading_structure("📊 Platform Overview", level=1, id="overview-title")
    # Add cache refresh button
    col_header, col_refresh = st.columns([10, 1])
    with col_header:
        st.header("📊 Overview")
    with col_refresh:
        if st.button("🔄", help="Refresh data", key="refresh_overview"):
            clear_all_caches()
            st.rerun()

    # Use cached data for metrics
    programs_data = fetch_programs()
    datasets_data = fetch_datasets(limit=1000)  # Higher limit for overview
    experiments_data = fetch_experiments(limit=1000)
    compounds_data = fetch_compounds(limit=1000)

    with db_session() as db:
        # Metrics row at top (using cached data where possible)
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            exp_count = len(experiments_data)  # Use cached count
            st.metric("Experiments", exp_count, delta=None)
        with col2:
            comp_count = len(compounds_data)  # Use cached count
            st.metric("Compounds", comp_count, delta=None)
        with col3:
            # Sample count still needs database query (not in cached data)
            sample_count = get_sample_count(db)
            st.metric("Samples", sample_count, delta=None)
        with col4:
            # Discovery count still needs database query (time-based)
            disc_count = get_discovery_count(days=7, db=db)
            st.metric("Discoveries (7d)", disc_count, delta=None)

        st.markdown("---")

        # Activity stats
        stats = get_activity_stats(db)

        # Stats row
        st.subheader("📈 Platform Statistics")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Experiments", stats["total_experiments"])
        with col2:
            st.metric("Total Compounds", stats["total_compounds"])
        with col3:
            st.metric("Total Datasets", stats["total_datasets"])
        with col4:
            st.metric("Total Users", stats["total_users"])

        st.markdown("---")

        # Quick Actions
        st.subheader("⚡ Quick Actions")
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("➕ New Experiment", use_container_width=True, type="primary"):
                st.session_state["selected_page"] = "Experiments"
                st.rerun()
        with col2:
            if st.button("⚗️ Register Compound", use_container_width=True, type="primary"):
                st.session_state["selected_page"] = "Chemistry"
                st.rerun()
        with col3:
            if st.button("🔬 Run Analysis", use_container_width=True, type="primary"):
                st.session_state["selected_page"] = "Analysis Tools"
                st.rerun()

        st.markdown("---")

        # Recent Activity
        st.subheader("🕐 Recent Activity")

        # Recent experiments
        recent_experiments = get_recent_experiments(db, limit=5)
        if recent_experiments:
            st.markdown("**Recent Experiments**")
            for exp in recent_experiments:
                created_by = f" by {exp['created_by']}" if exp['created_by'] else ""
                created_at = exp['created_at'][:10] if exp['created_at'] else "Unknown"
                # Create a link-like display (clicking would navigate to Experiments page)
                st.markdown(f"- **[{exp['name']}](?page=Experiments)** ({created_at}){created_by}")

        # Recent discoveries
        recent_discoveries = get_recent_discoveries(db, limit=5)
        if recent_discoveries:
            st.markdown("**Recent Discoveries**")
            for disc in recent_discoveries:
                st.markdown(f"- **{disc['study_id']}** ({disc['repository']}) - {disc['title'] or 'No title'} ({disc['discovered_at'][:10] if disc['discovered_at'] else 'Unknown'})")

        # Recent compounds
        recent_compounds = get_recent_compounds(db, limit=5)
        if recent_compounds:
            st.markdown("**Recent Compounds**")
            for comp in recent_compounds:
                st.markdown(f"- **{comp['compound_id']}** - {comp['smiles']} ({comp['created_at'][:10] if comp['created_at'] else 'Unknown'})")

        st.markdown("---")

        # Original statistics section
        # Statistics - Core entities
        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            dataset_count = len(datasets_data)  # Use cached count
            st.metric("Datasets", dataset_count)

        with col2:
            program_count = len(programs_data)  # Use cached count
            st.metric("Programs", program_count)

        with col3:
            experiment_count = len(experiments_data)  # Use cached count
            st.metric("Experiments", experiment_count)

        with col4:
            # Features still need database query (not in cached data)
            feature_count = db.query(func.count(Feature.id)).scalar()
            st.metric("Features", feature_count)

        with col5:
            # Signatures still need database query (not in cached data)
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
            compound_count = len(compounds_data)  # Use cached count
            st.metric("Compounds", compound_count)

        with col5:
            campaign_count = db.query(func.count(HTSCampaign.id)).scalar()
            st.metric("HTS Campaigns", campaign_count)

        # Datasets by Omics Type (using cached data)
        st.subheader("📈 Datasets by Omics Type")
        
        # Calculate omics type distribution from cached data
        omics_count_dict = {}
        for dataset in datasets_data:
            omics_type = dataset.get("omics_type", "Unknown")
            omics_count_dict[omics_type] = omics_count_dict.get(omics_type, 0) + 1
        
        omics_counts = list(omics_count_dict.items())

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
