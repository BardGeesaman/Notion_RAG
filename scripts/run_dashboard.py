#!/usr/bin/env python3
"""
Streamlit dashboard for the Amprenta Multi-Omics Platform.

This dashboard provides a visual interface for browsing and exploring
data stored in Postgres.

Usage:
    python scripts/run_dashboard.py
    # Or: streamlit run scripts/run_dashboard.py
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import streamlit as st
from sqlalchemy.orm import Session
from sqlalchemy import func, select

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import (
    Dataset,
    Program,
    Experiment,
    Feature,
    Signature,
)

# Page configuration
st.set_page_config(
    page_title="Amprenta Multi-Omics Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Title
st.title("🧬 Amprenta Multi-Omics Platform")
st.markdown("**Data Dashboard** - Browse and explore your multi-omics data")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Select Page",
    ["Overview", "Datasets", "Programs", "Experiments", "Features", "Signatures"],
)

# Helper function to get database session
@st.cache_resource
def get_database_session():
    """Get database session (cached)."""
    return next(get_db())


# Overview Page
if page == "Overview":
    st.header("📊 Overview")
    
    db = get_database_session()
    
    # Statistics
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
    
    # Datasets by Omics Type
    st.subheader("Datasets by Omics Type")
    omics_counts = (
        db.query(Dataset.omics_type, func.count(Dataset.id))
        .group_by(Dataset.omics_type)
        .all()
    )
    
    if omics_counts:
        import pandas as pd
        df_omics = pd.DataFrame(omics_counts, columns=["Omics Type", "Count"])
        
        col1, col2 = st.columns(2)
        with col1:
            st.bar_chart(df_omics.set_index("Omics Type"))
        with col2:
            st.dataframe(df_omics, use_container_width=True, hide_index=True)
            
            # Export button
            csv = df_omics.to_csv(index=False)
            st.download_button(
                label="📥 Download CSV",
                data=csv,
                file_name="omics_distribution.csv",
                mime="text/csv",
            )
    else:
        st.info("No datasets yet. Ingest some data to see statistics!")
        st.markdown("""
        **To get started:**
        1. Ingest a dataset: `python scripts/ingest_lipidomics.py --file data.csv --create-page`
        2. Refresh this page to see your data
        """)
    
    # Recent Datasets
    st.subheader("Recent Datasets")
    recent_datasets = (
        db.query(Dataset)
        .order_by(Dataset.created_at.desc())
        .limit(10)
        .all()
    )
    
    if recent_datasets:
        dataset_data = []
        for ds in recent_datasets:
            dataset_data.append({
                "Name": ds.name,
                "Omics Type": ds.omics_type,
                "Created": ds.created_at.strftime("%Y-%m-%d %H:%M"),
            })
        df_recent = pd.DataFrame(dataset_data)
        st.dataframe(df_recent, use_container_width=True)
    else:
        st.info("No datasets found.")


# Datasets Page
elif page == "Datasets":
    st.header("📁 Datasets")
    
    db = get_database_session()
    
    # Filters
    col1, col2 = st.columns(2)
    with col1:
        # Get distinct omics types
        omics_types = [row[0] for row in db.query(Dataset.omics_type).distinct().all() if row[0]]
        omics_filter = st.selectbox(
            "Filter by Omics Type",
            ["All"] + omics_types,
        )
    with col2:
        search_term = st.text_input("Search by name", "")
    
    # Query datasets
    query = db.query(Dataset)
    
    if omics_filter != "All":
        query = query.filter(Dataset.omics_type == omics_filter)
    
    if search_term:
        query = query.filter(Dataset.name.ilike(f"%{search_term}%"))
    
    datasets = query.order_by(Dataset.created_at.desc()).all()
    
    st.metric("Total Datasets", len(datasets))
    
    if datasets:
        # Summary table
        import pandas as pd
        summary_data = []
        for ds in datasets:
            summary_data.append({
                "Name": ds.name,
                "Omics Type": ds.omics_type,
                "Created": ds.created_at.strftime("%Y-%m-%d"),
                "Files": len(ds.file_paths) if ds.file_paths else 0,
            })
        df_summary = pd.DataFrame(summary_data)
        st.dataframe(df_summary, use_container_width=True, hide_index=True)
        
        # Export buttons
        col1, col2 = st.columns(2)
        with col1:
            csv_summary = df_summary.to_csv(index=False)
            st.download_button(
                label="📥 Download Summary (CSV)",
                data=csv_summary,
                file_name="datasets_summary.csv",
                mime="text/csv",
            )
        with col2:
            # Full dataset export
            full_data = []
            for ds in datasets:
                full_data.append({
                    "ID": str(ds.id),
                    "Name": ds.name,
                    "Omics Type": ds.omics_type,
                    "Description": ds.description or "",
                    "Created": ds.created_at.strftime("%Y-%m-%d %H:%M:%S"),
                    "Updated": ds.updated_at.strftime("%Y-%m-%d %H:%M:%S"),
                    "File Paths": "; ".join(ds.file_paths) if ds.file_paths else "",
                    "Disease": ", ".join(ds.disease) if ds.disease else "",
                    "Notion ID": ds.notion_page_id or "",
                })
            df_full = pd.DataFrame(full_data)
            csv_full = df_full.to_csv(index=False)
            st.download_button(
                label="📥 Download Full Data (CSV)",
                data=csv_full,
                file_name="datasets_full.csv",
                mime="text/csv",
            )
        
        st.markdown("---")
        st.subheader("Dataset Details")
        
        # Display datasets
        for dataset in datasets:
            with st.expander(f"**{dataset.name}** ({dataset.omics_type})"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{dataset.id}`")
                    st.write(f"**Omics Type:** {dataset.omics_type}")
                    if dataset.description:
                        st.write(f"**Description:** {dataset.description}")
                with col2:
                    st.write(f"**Created:** {dataset.created_at.strftime('%Y-%m-%d %H:%M')}")
                    st.write(f"**Updated:** {dataset.updated_at.strftime('%Y-%m-%d %H:%M')}")
                    if dataset.notion_page_id:
                        st.write(f"**Notion ID:** `{dataset.notion_page_id}`")
                
                if dataset.file_paths:
                    st.write("**File Paths:**")
                    for fp in dataset.file_paths:
                        st.code(fp, language=None)
                
                if dataset.disease:
                    st.write(f"**Disease:** {', '.join(dataset.disease)}")
                
                # Related entities (Postgres relationships)
                if dataset.programs:
                    st.write(f"**Linked Programs:** {len(dataset.programs)}")
                    program_names = [p.name for p in dataset.programs[:5]]
                    if program_names:
                        st.caption(f"Programs: {', '.join(program_names)}" + 
                                  (f" ... and {len(dataset.programs) - 5} more" if len(dataset.programs) > 5 else ""))
                
                if dataset.experiments:
                    st.write(f"**Linked Experiments:** {len(dataset.experiments)}")
                    experiment_names = [e.name for e in dataset.experiments[:5]]
                    if experiment_names:
                        st.caption(f"Experiments: {', '.join(experiment_names)}" + 
                                  (f" ... and {len(dataset.experiments) - 5} more" if len(dataset.experiments) > 5 else ""))
                
                # Linked features (Postgres)
                if dataset.features:
                    st.write(f"**Linked Features:** {len(dataset.features)}")
                    # Show feature types breakdown
                    feature_types = {}
                    for feature in dataset.features:
                        feature_types[feature.feature_type] = feature_types.get(feature.feature_type, 0) + 1
                    if feature_types:
                        type_breakdown = ", ".join([f"{ftype}: {count}" for ftype, count in feature_types.items()])
                        st.caption(f"Feature types: {type_breakdown}")
                    
                    # Show sample features (first 10)
                    if len(dataset.features) > 0:
                        sample_features = dataset.features[:10]
                        feature_names = [f.name for f in sample_features]
                        st.caption(f"Sample features: {', '.join(feature_names[:5])}" + 
                                  (f" ... and {len(dataset.features) - 5} more" if len(dataset.features) > 5 else ""))
    else:
        st.info("No datasets found matching your filters.")


# Programs Page
elif page == "Programs":
    st.header("🏥 Programs")
    
    db = get_database_session()
    
    search_term = st.text_input("Search programs by name", "")
    
    query = db.query(Program)
    if search_term:
        query = query.filter(Program.name.ilike(f"%{search_term}%"))
    
    programs = query.order_by(Program.created_at.desc()).all()
    
    st.metric("Total Programs", len(programs))
    
    if programs:
        for program in programs:
            with st.expander(f"**{program.name}**"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{program.id}`")
                    if program.description:
                        st.write(f"**Description:** {program.description}")
                with col2:
                    st.write(f"**Created:** {program.created_at.strftime('%Y-%m-%d %H:%M')}")
                    if program.disease:
                        st.write(f"**Disease:** {', '.join(program.disease)}")
                
                # Count related datasets
                dataset_count = len(program.datasets)
                if dataset_count > 0:
                    st.write(f"**Related Datasets:** {dataset_count}")
    else:
        st.info("No programs found.")


# Experiments Page
elif page == "Experiments":
    st.header("🔬 Experiments")
    
    db = get_database_session()
    
    search_term = st.text_input("Search experiments by name", "")
    
    query = db.query(Experiment)
    if search_term:
        query = query.filter(Experiment.name.ilike(f"%{search_term}%"))
    
    experiments = query.order_by(Experiment.created_at.desc()).all()
    
    st.metric("Total Experiments", len(experiments))
    
    if experiments:
        # Summary table
        import pandas as pd
        experiment_data = []
        for exp in experiments:
            experiment_data.append({
                "Name": exp.name,
                "Type": exp.type or "",
                "Disease": ", ".join(exp.disease) if exp.disease else "",
                "Matrix": ", ".join(exp.matrix) if exp.matrix else "",
                "Created": exp.created_at.strftime("%Y-%m-%d"),
                "Datasets": len(exp.datasets),
            })
        df_experiments = pd.DataFrame(experiment_data)
        st.dataframe(df_experiments, use_container_width=True, hide_index=True)
        
        # Export button
        csv_experiments = df_experiments.to_csv(index=False)
        st.download_button(
            label="📥 Download Experiments (CSV)",
            data=csv_experiments,
            file_name="experiments.csv",
            mime="text/csv",
        )
        
        st.markdown("---")
        st.subheader("Experiment Details")
        
        for experiment in experiments:
            with st.expander(f"**{experiment.name}**"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{experiment.id}`")
                    if experiment.type:
                        st.write(f"**Type:** {experiment.type}")
                    if experiment.description:
                        st.write(f"**Description:** {experiment.description}")
                with col2:
                    st.write(f"**Created:** {experiment.created_at.strftime('%Y-%m-%d %H:%M')}")
                    if experiment.disease:
                        st.write(f"**Disease:** {', '.join(experiment.disease)}")
                    if experiment.matrix:
                        st.write(f"**Matrix:** {', '.join(experiment.matrix)}")
                
                # Count related datasets
                dataset_count = len(experiment.datasets)
                if dataset_count > 0:
                    st.write(f"**Related Datasets:** {dataset_count}")
    else:
        st.info("No experiments found.")


# Features Page
elif page == "Features":
    st.header("🧪 Features")
    
    db = get_database_session()
    
    col1, col2 = st.columns(2)
    with col1:
        feature_type_filter = st.selectbox(
            "Filter by Feature Type",
            ["All", "gene", "protein", "metabolite", "lipid"],
        )
    with col2:
        search_term = st.text_input("Search features by name", "")
    
    query = db.query(Feature)
    
    if feature_type_filter != "All":
        query = query.filter(Feature.feature_type == feature_type_filter)
    
    if search_term:
        query = query.filter(Feature.name.ilike(f"%{search_term}%"))
    
    features = query.order_by(Feature.name).limit(100).all()
    
    st.metric("Features Found", len(features))
    if len(features) == 100:
        st.info("Showing first 100 features. Use filters to narrow down.")
    
    if features:
        # Group by feature type
        import pandas as pd
        feature_data = []
        for feature in features:
            # Count datasets linked to this feature
            dataset_count = len(feature.datasets) if hasattr(feature, 'datasets') and feature.datasets else 0
            
            feature_data.append({
                "Name": feature.name,
                "Type": feature.feature_type,
                "Normalized Name": feature.normalized_name or "-",
                "Linked Datasets": dataset_count,
            })
        df_features = pd.DataFrame(feature_data)
        st.dataframe(df_features, use_container_width=True)
        
        # Feature type distribution
        if len(df_features) > 0:
            type_counts = df_features["Type"].value_counts()
            if len(type_counts) > 1:
                st.subheader("Feature Type Distribution")
                col1, col2 = st.columns(2)
                with col1:
                    st.bar_chart(type_counts)
                with col2:
                    st.dataframe(type_counts.reset_index().rename(columns={"index": "Type", "Type": "Count"}))
        
        # Export button
        csv_features = df_features.to_csv(index=False)
        st.download_button(
            label="📥 Download Features (CSV)",
            data=csv_features,
            file_name="features.csv",
            mime="text/csv",
        )
    else:
        st.info("No features found.")


# Signatures Page
elif page == "Signatures":
    st.header("📋 Signatures")
    
    db = get_database_session()
    
    search_term = st.text_input("Search signatures by name", "")
    
    query = db.query(Signature)
    if search_term:
        query = query.filter(Signature.name.ilike(f"%{search_term}%"))
    
    signatures = query.order_by(Signature.created_at.desc()).all()
    
    st.metric("Total Signatures", len(signatures))
    
    if signatures:
        # Summary table
        import pandas as pd
        signature_data = []
        for sig in signatures:
            signature_data.append({
                "Name": sig.name,
                "Description": sig.description or "",
                "Modalities": ", ".join(sig.modalities) if sig.modalities else "",
                "Components": len(sig.components),
                "Created": sig.created_at.strftime("%Y-%m-%d"),
            })
        df_signatures = pd.DataFrame(signature_data)
        st.dataframe(df_signatures, use_container_width=True, hide_index=True)
        
        # Modalities distribution
        if len(df_signatures) > 0 and any(df_signatures["Modalities"]):
            st.subheader("Signature Modalities Distribution")
            # Count modalities
            all_modalities = []
            for mods in df_signatures["Modalities"]:
                if mods:
                    all_modalities.extend([m.strip() for m in mods.split(",")])
            if all_modalities:
                modality_counts = pd.Series(all_modalities).value_counts()
                st.bar_chart(modality_counts)
        
        # Export button
        csv_signatures = df_signatures.to_csv(index=False)
        st.download_button(
            label="📥 Download Signatures (CSV)",
            data=csv_signatures,
            file_name="signatures.csv",
            mime="text/csv",
        )
        
        st.markdown("---")
        st.subheader("Signature Details")
        
        for signature in signatures:
            with st.expander(f"**{signature.name}**"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{signature.id}`")
                    if signature.description:
                        st.write(f"**Description:** {signature.description}")
                    if signature.modalities:
                        st.write(f"**Modalities:** {', '.join(signature.modalities)}")
                with col2:
                    st.write(f"**Created:** {signature.created_at.strftime('%Y-%m-%d %H:%M')}")
                    component_count = len(signature.components)
                    st.write(f"**Components:** {component_count}")
    else:
        st.info("No signatures found.")


# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("**Amprenta Multi-Omics Platform**")
st.sidebar.markdown("Data stored in Postgres")

# Add refresh button
if st.sidebar.button("🔄 Refresh Data"):
    st.cache_resource.clear()
    st.rerun()

# Add API link
st.sidebar.markdown("---")
st.sidebar.markdown("**API Access**")
st.sidebar.markdown("[FastAPI Docs](http://localhost:8000/docs)")
st.sidebar.markdown("[API Health](http://localhost:8000/health)")

# Add quick actions
st.sidebar.markdown("---")
st.sidebar.markdown("**Quick Actions**")
st.sidebar.markdown("""
- Ingest data via CLI scripts
- Use FastAPI for programmatic access
- Query Postgres directly for advanced queries
""")

