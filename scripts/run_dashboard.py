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
import os

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import streamlit as st
from uuid import UUID
from amprenta_rag.auth.session import (
    get_current_user,
    clear_session,
    is_authenticated,
    update_last_activity,
    check_session_timeout,
    get_session_remaining,
)
from amprenta_rag.auth.audit import log_logout
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import UserFavorite, Experiment, Compound, Signature
from scripts.dashboard.help_content import get_help, search_help
from scripts.dashboard.help_chat import render_help_chat
from scripts.dashboard.themes import apply_theme, get_theme_names
from scripts.dashboard.shortcuts import inject_shortcuts_js, handle_shortcut_action, render_shortcuts_help
from amprenta_rag.notifications.service import (
    get_unread_count,
    get_unread_notifications,
    mark_as_read,
    mark_all_read,
)
from amprenta_rag.utils.global_search import global_search
from amprenta_rag.utils.bookmarks import get_user_bookmarks
AUTH_DISABLED = os.environ.get("DISABLE_AUTH", "").lower() in ("1", "true", "yes")

# Page groups
DISCOVERY_PAGES = ["Overview", "Experiments", "Discovery Workflow"]
ANALYSIS_PAGES = ["Analysis Tools", "Chemistry", "Visualizations"]
ELN_PAGES = ["Protocols", "Sample Inventory", "Q&A Tracker"]
ADMIN_PAGES = ["Audit Logs", "Teams & Projects", "Feedback", "Import Data"]
ALL_PAGES = [
    "Overview",
    "Getting Started",
    "Evaluation Wizard",
    "Chat",
    "Lab Notebook",
    "Sample Inventory",
    "Search",
    "Data Ingestion",
    "Repositories",
    "Discovery Workflow",
    "Analysis Tools",
    "Visualizations",
    "Quality Checks",
    "Statistical Analysis",
    "Discovery",
    "Coverage Map",
    "Feature Recurrence",
    "Evidence Report",
    "Data Management",
    "System Health",
    "Audit Logs",
    "Relationships",
    "Datasets",
    "Programs",
    "Experiments",
    "Features",
    "Signatures",
    "Protocols",
    "Q&A Tracker",
    "Teams & Projects",
    "Feedback",
    "Literature",
    "Emails",
    "RAG Chunks",
    "Chemistry",
    "RAG Query",
    "Cross-Omics",
    "Import Data",
    "Compare",
    "Timeline",
    "Data Quality",
]

# LAZY IMPORTS: Don't import pages until needed to avoid cascading import failures
# from scripts.dashboard.pages import ...

# Page configuration
st.set_page_config(
    page_title="Amprenta Multi-Omics Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Update last activity timestamp
if is_authenticated():
    update_last_activity()

# Authentication check
if not AUTH_DISABLED and not is_authenticated():
    from scripts.dashboard.pages.auth.login import render_login_page

    render_login_page()
    st.stop()

# Check session timeout
if not AUTH_DISABLED and is_authenticated():
    if check_session_timeout(timeout_minutes=30):
        clear_session()
        st.warning("⏱️ Your session has expired. Please log in again.")
        from scripts.dashboard.pages.auth.login import render_login_page
        render_login_page()
        st.stop()
    else:
        # Show warning if less than 5 minutes remaining
        remaining = get_session_remaining(timeout_minutes=30)
        if 0 < remaining < 5:
            st.warning(f"⏱️ Your session will expire in {remaining} minute{'s' if remaining != 1 else ''}. Please save your work.")

# Mock user when auth disabled (e.g., testing)
if AUTH_DISABLED and not is_authenticated():
    from amprenta_rag.auth.session import set_current_user

    set_current_user({"id": "test", "username": "dev", "email": "dev@local", "role": "admin"})

# Admin registration page
if st.session_state.get("show_register"):
    from scripts.dashboard.pages.auth.register import render_register_page

    render_register_page()
    if st.button("← Back to Dashboard"):
        st.session_state.pop("show_register", None)
        st.rerun()
    st.stop()

# Title
st.title("🧬 Amprenta Multi-Omics Platform")
st.markdown("**Data Dashboard** - Browse and explore your multi-omics data")

# Sidebar navigation
user = get_current_user()

# Initialize theme
if "theme" not in st.session_state:
    st.session_state["theme"] = "dark"

# Apply theme (must be before sidebar content)
apply_theme(st.session_state["theme"])

if user:
    with st.sidebar:
        st.markdown(f"**👤 {user['username']}** ({user['role']})")
        if st.button("Logout", key="logout_btn"):
            log_logout(user.get("id"), user.get("username"))
            clear_session()
            st.rerun()
        st.divider()
        
        # Notifications
        user_id = user.get("id")
        if user_id and user_id != "test":
            try:
                db_gen = get_db()
                db = next(db_gen)
                try:
                    unread_count = get_unread_count(user_id, db)
                    notifications = get_unread_notifications(user_id, db, limit=10)
                    
                    # Notification bell with count
                    if unread_count > 0:
                        st.markdown(f"**🔔 Notifications ({unread_count})**")
                    else:
                        st.markdown("**🔔 Notifications**")
                    
                    # Notifications expander
                    with st.expander("View Notifications", expanded=False):
                        if notifications:
                            for notif in notifications:
                                col1, col2 = st.columns([3, 1])
                                with col1:
                                    st.markdown(f"**{notif.title}**")
                                    if notif.message:
                                        st.caption(notif.message)
                                    st.caption(f"{notif.created_at.strftime('%Y-%m-%d %H:%M') if notif.created_at else ''} • {notif.notification_type}")
                                with col2:
                                    if st.button("✓", key=f"mark_read_{notif.id}", help="Mark as read"):
                                        mark_as_read(str(notif.id), db)
                                        st.rerun()
                                st.divider()
                            
                            if st.button("Mark All Read", key="mark_all_read", use_container_width=True):
                                mark_all_read(user_id, db)
                                st.rerun()
                        else:
                            st.info("No unread notifications")
                finally:
                    db_gen.close()
            except Exception as e:
                st.error(f"Error loading notifications: {e}")
        
        st.divider()
        
        # Theme selector
        theme_options = get_theme_names()
        selected_theme = st.selectbox(
            "Theme",
            theme_options,
            index=theme_options.index(st.session_state["theme"]) if st.session_state["theme"] in theme_options else 0,
            key="theme_selector"
        )
        if selected_theme != st.session_state["theme"]:
            st.session_state["theme"] = selected_theme
            st.rerun()
        
        st.divider()
        if user.get("role") == "admin":
            if st.button("➕ Register User", key="register_btn"):
                from scripts.dashboard.pages.auth.register import render_register_page

                st.session_state["show_register"] = True
                st.rerun()
else:
    # Theme selector for non-authenticated users
    with st.sidebar:
        theme_options = get_theme_names()
        selected_theme = st.selectbox(
            "Theme",
            theme_options,
            index=theme_options.index(st.session_state["theme"]) if st.session_state["theme"] in theme_options else 0,
            key="theme_selector"
        )
        if selected_theme != st.session_state["theme"]:
            st.session_state["theme"] = selected_theme
            st.rerun()

# Initialize recent pages
if "recent_pages" not in st.session_state:
    st.session_state["recent_pages"] = []

# Helper functions for favorites
def get_user_favorites(user_id: str) -> list[str]:
    """Get user's favorite pages."""
    if not user_id or user_id == "test":
        return []
    try:
        db_gen = get_db()
        db = next(db_gen)
        try:
            favorites = db.query(UserFavorite).filter(
                UserFavorite.user_id == UUID(user_id)
            ).order_by(UserFavorite.created_at.desc()).all()
            return [f.page_name for f in favorites]
        finally:
            db_gen.close()
    except Exception:
        return []

def toggle_favorite(user_id: str, page_name: str) -> None:
    """Add or remove favorite."""
    if not user_id or user_id == "test":
        return
    try:
        db_gen = get_db()
        db = next(db_gen)
        try:
            existing = db.query(UserFavorite).filter(
                UserFavorite.user_id == UUID(user_id),
                UserFavorite.page_name == page_name
            ).first()
            if existing:
                db.delete(existing)
            else:
                favorite = UserFavorite(
                    user_id=UUID(user_id),
                    page_name=page_name
                )
                db.add(favorite)
            db.commit()
        finally:
            db_gen.close()
    except Exception:
        pass

def update_recent_pages(page_name: str) -> None:
    """Update recent pages list."""
    if page_name not in st.session_state["recent_pages"]:
        st.session_state["recent_pages"].insert(0, page_name)
        st.session_state["recent_pages"] = st.session_state["recent_pages"][:5]
    else:
        # Move to front
        st.session_state["recent_pages"].remove(page_name)
        st.session_state["recent_pages"].insert(0, page_name)

# Sidebar navigation
st.sidebar.title("Navigation")

# Global search
search_query = st.sidebar.text_input("🔍 Search...", key="global_search", placeholder="Search experiments, compounds, datasets...")
if search_query:
    with st.sidebar.expander("🔍 Search Results", expanded=True):
        try:
            db_gen = get_db()
            db = next(db_gen)
            try:
                results = global_search(search_query, db, limit=5)
                
                # Display results grouped by type
                if results["experiments"]:
                    st.markdown("**Experiments**")
                    for exp in results["experiments"]:
                        if st.button(f"🔬 {exp['name']}", key=f"search_exp_{exp['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Experiments"
                            st.session_state["selected_experiment_id"] = exp["id"]
                            st.rerun()
                
                if results["compounds"]:
                    st.markdown("**Compounds**")
                    for comp in results["compounds"]:
                        if st.button(f"⚗️ {comp['compound_id']}", key=f"search_comp_{comp['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Chemistry"
                            st.session_state["selected_compound_id"] = comp["id"]
                            st.rerun()
                
                if results["signatures"]:
                    st.markdown("**Signatures**")
                    for sig in results["signatures"]:
                        if st.button(f"📊 {sig['name']}", key=f"search_sig_{sig['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Signatures"
                            st.session_state["selected_signature_id"] = sig["id"]
                            st.rerun()
                
                if results["datasets"]:
                    st.markdown("**Datasets**")
                    for ds in results["datasets"]:
                        if st.button(f"📁 {ds['name']}", key=f"search_ds_{ds['id']}", use_container_width=True):
                            st.session_state["selected_page"] = "Datasets"
                            st.session_state["selected_dataset_id"] = ds["id"]
                            st.rerun()
                
                if not any([results["experiments"], results["compounds"], results["signatures"], results["datasets"]]):
                    st.info("No results found")
            finally:
                db_gen.close()
        except Exception as e:
            st.error(f"Search error: {e}")

st.sidebar.divider()
page = None

# Favorites section
if user and user.get("id") != "test":
    favorites = get_user_favorites(user.get("id"))
    if favorites:
        with st.sidebar.expander("⭐ Favorites", expanded=True):
            for fav_page in favorites:
                if st.button(f"📌 {fav_page}", key=f"fav_{fav_page}", use_container_width=True):
                    page = fav_page
                    st.session_state["selected_page"] = page
                    update_recent_pages(page)
                    st.rerun()
        st.sidebar.divider()

# Bookmarks section
if user and user.get("id") != "test":
    try:
        db_gen = get_db()
        db = next(db_gen)
        try:
            bookmarks = get_user_bookmarks(user.get("id"), db)
            if bookmarks:
                with st.sidebar.expander("📌 Bookmarks", expanded=False):
                    # Group by entity type
                    experiments = [b for b in bookmarks if b.entity_type == "experiment"]
                    compounds = [b for b in bookmarks if b.entity_type == "compound"]
                    signatures = [b for b in bookmarks if b.entity_type == "signature"]
                    
                    if experiments:
                        st.markdown("**Experiments**")
                        for bm in experiments:
                            exp = db.query(Experiment).filter(Experiment.id == bm.entity_id).first()
                            if exp:
                                if st.button(f"🔬 {exp.name}", key=f"bm_exp_{bm.id}", use_container_width=True):
                                    st.session_state["selected_page"] = "Experiments"
                                    st.session_state["selected_experiment_id"] = str(bm.entity_id)
                                    st.rerun()
                    
                    if compounds:
                        st.markdown("**Compounds**")
                        for bm in compounds:
                            comp = db.query(Compound).filter(Compound.id == bm.entity_id).first()
                            if comp:
                                if st.button(f"⚗️ {comp.compound_id}", key=f"bm_comp_{bm.id}", use_container_width=True):
                                    st.session_state["selected_page"] = "Chemistry"
                                    st.session_state["selected_compound_id"] = str(bm.entity_id)
                                    st.rerun()
                    
                    if signatures:
                        st.markdown("**Signatures**")
                        for bm in signatures:
                            sig = db.query(Signature).filter(Signature.id == bm.entity_id).first()
                            if sig:
                                if st.button(f"📊 {sig.name}", key=f"bm_sig_{bm.id}", use_container_width=True):
                                    st.session_state["selected_page"] = "Signatures"
                                    st.session_state["selected_signature_id"] = str(bm.entity_id)
                                    st.rerun()
        finally:
            db_gen.close()
    except Exception:
        pass  # Silently fail if bookmarks can't be loaded

# Recent pages
if st.session_state.get("recent_pages"):
    with st.sidebar.expander("🕐 Recent", expanded=False):
        for recent_page in st.session_state["recent_pages"][:5]:
            if st.button(f"🕐 {recent_page}", key=f"recent_{recent_page}", use_container_width=True):
                page = recent_page
                st.session_state["selected_page"] = page
                update_recent_pages(page)
                st.rerun()
    st.sidebar.divider()

# Initialize selected page from session state
if "selected_page" not in st.session_state:
    st.session_state["selected_page"] = "Overview"

# Handle shortcut actions (check for navigation triggers)
shortcut_actions = [
    "navigate_home", "navigate_experiments", "navigate_chemistry", 
    "navigate_analysis", "navigate_search", "navigate_ingestion",
    "navigate_protocols", "navigate_qa"
]
for action in shortcut_actions:
    if st.session_state.get(f"shortcut_{action}", False):
        handle_shortcut_action(action)
        st.session_state.pop(f"shortcut_{action}", None)
        break

# Show shortcuts help if requested
if st.session_state.get("show_shortcuts", False):
    render_shortcuts_help()
    if st.button("Close", key="close_shortcuts"):
        st.session_state.pop("show_shortcuts", None)
        st.rerun()

# Command palette
from scripts.dashboard.components.command_palette import render_command_palette

render_command_palette()

# Command palette
from scripts.dashboard.components.command_palette import render_command_palette
render_command_palette()

# Page groups
if page is None:
    page = st.session_state.get("selected_page", "Overview")
    
    with st.sidebar.expander("🔍 Discovery", expanded=True):
        for p in DISCOVERY_PAGES:
            if st.button(p, key=f"btn_{p}", use_container_width=True):
                page = p
                st.session_state["selected_page"] = p
                update_recent_pages(p)
                st.rerun()
    
    with st.sidebar.expander("📊 Analysis", expanded=False):
        for p in ANALYSIS_PAGES:
            if st.button(p, key=f"btn_{p}", use_container_width=True):
                page = p
                st.session_state["selected_page"] = p
                update_recent_pages(p)
                st.rerun()
    
    with st.sidebar.expander("📋 ELN", expanded=False):
        for p in ELN_PAGES:
            if st.button(p, key=f"btn_{p}", use_container_width=True):
                page = p
                st.session_state["selected_page"] = p
                update_recent_pages(p)
                st.rerun()
    
    if user and user.get("role") == "admin":
        with st.sidebar.expander("⚙️ Admin", expanded=False):
            for p in ADMIN_PAGES:
                if st.button(p, key=f"btn_{p}", use_container_width=True):
                    page = p
                    st.session_state["selected_page"] = p
                    update_recent_pages(p)
                    st.rerun()
    
    # Other pages
    other_pages = [p for p in ALL_PAGES if p not in DISCOVERY_PAGES + ANALYSIS_PAGES + ELN_PAGES + ADMIN_PAGES]
    with st.sidebar.expander("📚 Other Pages", expanded=False):
        for p in other_pages:
            if st.button(p, key=f"btn_{p}", use_container_width=True):
                page = p
                st.session_state["selected_page"] = p
                update_recent_pages(p)
                st.rerun()

# Get final page selection from session state
page = st.session_state.get("selected_page", "Overview")

# Star button for favorites
if user and user.get("id") != "test":
    favorites = get_user_favorites(user.get("id"))
    is_favorite = page in favorites
    star_label = "⭐" if is_favorite else "☆"
    if st.sidebar.button(f"{star_label} Favorite", key="toggle_favorite", use_container_width=True):
        toggle_favorite(user.get("id"), page)
        st.rerun()

# Help section
with st.sidebar.expander("❓ Help", expanded=False):
    # Search box
    help_query = st.text_input("Search help", key="help_search")
    if help_query:
        results = search_help(help_query)
        for r in results[:3]:
            st.markdown(f"**{r['title']}**: {r['description']}")
    else:
        # Show help for current page
        current_help = get_help(page)
        if current_help:
            st.markdown(f"### {current_help['title']}")
            st.write(current_help['description'])
            if current_help.get('tips'):
                st.markdown("**Tips:**")
                for tip in current_help['tips']:
                    st.markdown(f"• {tip}")

# Help chat assistant
with st.sidebar.expander("🤖 Ask Assistant", expanded=False):
    render_help_chat()

# Update recent pages
update_recent_pages(page)

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
    elif page == "Sample Inventory":
        from scripts.dashboard.pages.sample_inventory import render_sample_inventory_page

        render_sample_inventory_page()
    elif page == "Search":
        from scripts.dashboard.pages.search import render_search_page

        render_search_page()
    elif page == "Data Ingestion":
        from scripts.dashboard.pages.ingestion import render_ingestion_page

        render_ingestion_page()
    elif page == "Repositories":
        from scripts.dashboard.pages.repositories import render_repositories_page

        render_repositories_page()
    elif page == "Discovery Workflow":
        from scripts.dashboard.pages.discovery_workflow import render_discovery_workflow_page

        render_discovery_workflow_page()
    elif page == "Analysis Tools":
        from scripts.dashboard.pages.analysis import render_analysis_page

        render_analysis_page()
    elif page == "Visualizations":
        from scripts.dashboard.pages.visualizations import render_visualizations_page

        render_visualizations_page()
    elif page == "Quality Checks":
        from scripts.dashboard.pages.quality_checks import render_quality_checks_page

        render_quality_checks_page()
    elif page == "Statistical Analysis":
        from scripts.dashboard.pages.statistical_analysis import render_statistical_analysis_page

        render_statistical_analysis_page()
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
        from scripts.dashboard.pages.system_health import render_system_health_page

        render_system_health_page()
    elif page == "Audit Logs":
        from scripts.dashboard.pages.audit_logs import render_audit_logs_page

        render_audit_logs_page()
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
    elif page == "Protocols":
        from scripts.dashboard.pages.protocols import render_protocols_page

        render_protocols_page()
    elif page == "Q&A Tracker":
        from scripts.dashboard.pages.qa_tracker import render_qa_tracker_page

        render_qa_tracker_page()
    elif page == "Teams & Projects":
        from scripts.dashboard.pages.teams import render_teams_page

        render_teams_page()
    elif page == "Feedback":
        from scripts.dashboard.pages.feedback import render_feedback_page

        render_feedback_page()
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
    elif page == "Import Data":
        from scripts.dashboard.pages.import_data import render_import_page

        render_import_page()
    elif page == "Compare":
        from scripts.dashboard.pages.compare import render_compare_page

        render_compare_page()
    elif page == "Timeline":
        from scripts.dashboard.pages.timeline import render_timeline_page

        render_timeline_page()
    elif page == "Data Quality":
        from scripts.dashboard.pages.data_quality import render_data_quality_page

        render_data_quality_page()
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
st.sidebar.caption("Press `?` for keyboard shortcuts")
st.sidebar.caption("Press `?` for keyboard shortcuts")

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
