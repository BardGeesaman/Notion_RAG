import streamlit as st

from amprenta_rag.database.models import Dataset
from scripts.dashboard.db_session import db_session


def render_getting_started_page():
    st.title("🚀 Getting Started with Amprenta Platform")
    st.markdown(
        """
    Welcome! This brief guide will walk you through the core workflows on the platform.

    👈 **Use the sidebar navigation** to access different pages.
    """
    )
    st.markdown("---")
    st.header("Quick Start Steps")
    cols = st.columns(4)
    # Step 1: Upload
    with cols[0]:
        st.subheader("1. Upload a Dataset")
        st.write("Ingest raw omics data with rich metadata.")
        st.info("📤 **Go to:** Data Ingestion (in sidebar)")
    # Step 2: QC & Coverage
    with cols[1]:
        st.subheader("2. View Your Data")
        st.write("Browse datasets, features, and metadata.")
        st.info("🗂️ **Go to:** Datasets, Features, or Programs (in sidebar)")
    # Step 3: Discovery/Signatures
    with cols[2]:
        st.subheader("3. Discover Signatures")
        st.write("Analyze for significant signals with signature discovery.")
        st.info("🔍 **Go to:** Signatures (in sidebar)")
    # Step 4: Ask the AI
    with cols[3]:
        st.subheader("4. Ask the Assistant")
        st.write("Let the AI explain, summarize, or guide you.")
        st.info("💬 **Go to:** Chat (in sidebar)")
    st.markdown("---")

    # Show some quick stats
    st.header("📊 Platform Overview")
    with db_session() as db:
        dataset_count = db.query(Dataset).count()
        st.metric("Total Datasets", dataset_count)

        if dataset_count > 0:
            st.success("✅ Your platform has data! Explore it using the sidebar navigation.")
        else:
            st.warning("⚠️ No datasets found. Start by ingesting data via the **Data Ingestion** page.")

    st.markdown("---")
    st.header("🎯 Key Features")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(
            """
        **Data Management:**
        - Browse programs, experiments, and datasets
        - View features across omics types
        - Explore relationships between entities
        """
        )
    with col2:
        st.markdown(
            """
        **Analysis & Query:**
        - RAG-powered literature search
        - Cross-omics reasoning
        - Signature discovery and validation
        """
        )
