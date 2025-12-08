#!/usr/bin/env python3
"""
Amprenta Multi-Omics Platform - Streamlit Dashboard

This is the PRIMARY USER INTERFACE for the Amprenta platform. It provides a comprehensive
web-based interface for all platform functionality.

**Architecture**: Postgres-First, Dashboard-Centric
- All data stored in PostgreSQL (sole system of record)
- 19 pages covering ingestion, search, analysis, and lab notebook
- No Notion dependency for core functionality

**Dashboard Pages** (19 total):
1. Overview - System statistics and health
2. Getting Started - Quick start guide
3. Evaluation Wizard - Dataset evaluation workflow
4. Chat - AI-powered chat interface
5. Lab Notebook - Electronic Lab Notebook (ELN)
6. Search - Global entity search
7. Data Ingestion - Upload datasets, signatures, chemistry data
8. Repositories - Import from public repositories (MW, GEO, PRIDE, MetaboLights)
9. Analysis Tools - Pathway enrichment, signature scoring
10. Discovery - Automated signature discovery (experimental)
11. Coverage Map - Program × Signature coverage matrices
12. Feature Recurrence - Feature co-occurrence analysis
13. Evidence Report - Cross-omics evidence summaries
14. Data Management - Link entities, edit metadata, bulk operations
15. System Health - Database status, configuration checks
16. Relationships - Visualize entity relationships
17-21. Browse pages: Datasets, Programs, Experiments, Features, Signatures
22-24. Content pages: Literature, Emails, RAG Chunks
25. Chemistry - Compounds, HTS campaigns, biochemical results
26. RAG Query - Semantic search with filtering
27. Cross-Omics - Multi-omics reasoning

**Dashboard Features**:
- Postgres-backed (fast SQL queries)
- Real-time data updates
- File upload/download
- Interactive visualizations
- Export to CSV/JSON
- Entity linking and relationships
- Search and filtering

Usage:
    # Development
    python scripts/run_dashboard.py
    # Or: streamlit run scripts/run_dashboard.py
    # Opens at: http://localhost:8501
    
    # Production (see docs/DEPLOYMENT_GUIDE.md)
    # Use systemd service + nginx reverse proxy
    
Documentation:
    - User workflows: docs/USER_GUIDE.md
    - Deployment: docs/DEPLOYMENT_GUIDE.md
    - Development: docs/DEVELOPER_GUIDE.md

Note:
    - Dashboard pages use lazy imports to prevent cascading failures
    - Database sessions managed via scripts/dashboard/db_session.py
    - Uses module import pattern for SQLAlchemy models (see docs/DEVELOPER_GUIDE.md)
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import streamlit as st

# LAZY IMPORTS: Don't import pages until needed to avoid cascading import failures
# from scripts.dashboard.pages import ...

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
    [
        "Overview",
        "Getting Started",
        "Evaluation Wizard",
        "Chat",
        "Lab Notebook",
        "Search",
        "Data Ingestion",
        "Repositories",
        "Analysis Tools",
        "Visualizations",
        "Quality Checks",
        "Discovery",
        "Coverage Map",
        "Feature Recurrence",
        "Evidence Report",
        "Data Management",
        "System Health",
        "Relationships",
        "Datasets",
        "Programs",
        "Experiments",
        "Features",
        "Signatures",
        "Literature",
        "Emails",
        "RAG Chunks",
        "Chemistry",
        "RAG Query",
        "Cross-Omics",
    ],
)

# Route to appropriate page (LAZY IMPORTS to prevent cascading failures)
try:
    if page == "Overview":
        from scripts.dashboard.pages.overview import render_overview_page

        render_overview_page()
    elif page == "Getting Started":
        from scripts.dashboard.pages.getting_started import render_getting_started_page

        render_getting_started_page()
    elif page == "Evaluation Wizard":
        from scripts.dashboard.pages.evaluation_wizard import render_evaluation_wizard

        render_evaluation_wizard()
    elif page == "Chat":
        from scripts.dashboard.pages.chat import render_chat_page

        render_chat_page()
    elif page == "Lab Notebook":
        from scripts.dashboard.pages.lab_notebook import render_lab_notebook_page

        render_lab_notebook_page()
    elif page == "Search":
        from scripts.dashboard.pages.search import render_search_page

        render_search_page()
    elif page == "Data Ingestion":
        from scripts.dashboard.pages.ingestion import render_ingestion_page

        render_ingestion_page()
    elif page == "Repositories":
        from scripts.dashboard.pages.repositories import render_repositories_page

        render_repositories_page()
    elif page == "Analysis Tools":
        from scripts.dashboard.pages.analysis import render_analysis_page

        render_analysis_page()
    elif page == "Visualizations":
        from scripts.dashboard.pages.visualizations import render_visualizations_page

        render_visualizations_page()
    elif page == "Quality Checks":
        from scripts.dashboard.pages.quality_checks import render_quality_checks_page

        render_quality_checks_page()
    elif page == "Discovery":
        from scripts.dashboard.pages.discovery import render_discovery_page

        render_discovery_page()
    elif page == "Coverage Map":
        from scripts.dashboard.pages.coverage import render_coverage_page

        render_coverage_page()
    elif page == "Feature Recurrence":
        from scripts.dashboard.pages.feature_recurrence import render_feature_recurrence_page

        render_feature_recurrence_page()
    elif page == "Evidence Report":
        from scripts.dashboard.pages.evidence_report import render_evidence_report_page

        render_evidence_report_page()
    elif page == "Data Management":
        from scripts.dashboard.pages.management import render_management_page

        render_management_page()
    elif page == "System Health":
        from scripts.dashboard.pages.health import render_health_page

        render_health_page()
    elif page == "Relationships":
        from scripts.dashboard.pages.relationships import render_relationships_page

        render_relationships_page()
    elif page == "Datasets":
        from scripts.dashboard.pages.datasets import render_datasets_page

        render_datasets_page()
    elif page == "Programs":
        from scripts.dashboard.pages.programs import render_programs_page

        render_programs_page()
    elif page == "Experiments":
        from scripts.dashboard.pages.experiments import render_experiments_page

        render_experiments_page()
    elif page == "Features":
        from scripts.dashboard.pages.features import render_features_page

        render_features_page()
    elif page == "Signatures":
        from scripts.dashboard.pages.signatures import render_signatures_page

        render_signatures_page()
    elif page == "Literature":
        from scripts.dashboard.pages.literature import render_literature_page

        render_literature_page()
    elif page == "Emails":
        from scripts.dashboard.pages.emails import render_emails_page

        render_emails_page()
    elif page == "RAG Chunks":
        from scripts.dashboard.pages.rag_chunks import render_rag_chunks_page

        render_rag_chunks_page()
    elif page == "Chemistry":
        from scripts.dashboard.pages.chemistry import render_chemistry_page

        render_chemistry_page()
    elif page == "RAG Query":
        from scripts.dashboard.pages.rag_query import render_rag_query_page

        render_rag_query_page()
    elif page == "Cross-Omics":
        from scripts.dashboard.pages.cross_omics import render_cross_omics_page

        render_cross_omics_page()
except ImportError as e:
    st.error(f"❌ Error loading page: {page}")
    st.error(f"Import error: {str(e)}")
    st.info("This page may depend on modules that are currently unavailable. Try another page.")
    with st.expander("Show full error details"):
        st.exception(e)
except Exception as e:
    st.error(f"❌ Error rendering page: {page}")
    st.exception(e)


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
st.sidebar.markdown(
    """
- Ingest data via CLI scripts
- Use FastAPI for programmatic access
- Query Postgres directly for advanced queries
"""
)
