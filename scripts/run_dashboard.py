#!/usr/bin/env python3
"""Amprenta Multi-Omics Platform - Streamlit Dashboard (modular)."""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

# Ensure project root on path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.dashboard.core import (
    AUTH_DISABLED,
    ALL_PAGES,
    DISCOVERY_PAGES,
    ANALYSIS_PAGES,
    ELN_PAGES,
    ADMIN_PAGES,
    check_authentication,
    handle_session_timeout,
    render_sidebar,
    route_to_page,
)
from amprenta_rag.auth.feature_permissions import get_visible_pages
from amprenta_rag.database.session import db_session


def main() -> None:
    st.set_page_config(
        page_title="Amprenta Multi-Omics Platform",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    if st.session_state.get("show_register"):
        from scripts.dashboard.pages.auth.register import render_register_page

        render_register_page()
        if st.button("← Back to Dashboard"):
            st.session_state.pop("show_register", None)
            st.rerun()
        st.stop()

    user = check_authentication(AUTH_DISABLED)
    handle_session_timeout(AUTH_DISABLED, timeout_minutes=30)

    # Visible pages by role
    user_role = user.get("role", "viewer") if user else "viewer"
    try:
        with db_session() as db:
            visible_pages = get_visible_pages(user_role, db)
    except Exception:
        visible_pages = ALL_PAGES
    if not visible_pages:
        visible_pages = ALL_PAGES

    st.title("🧬 Amprenta Multi-Omics Platform")
    st.markdown("**Data Dashboard** - Browse and explore your multi-omics data")

    page = render_sidebar(
        user,
        visible_pages,
        groups=type(
            "Groups",
            (),
            {
                "DISCOVERY_PAGES": DISCOVERY_PAGES,
                "ANALYSIS_PAGES": ANALYSIS_PAGES,
                "ELN_PAGES": ELN_PAGES,
                "ADMIN_PAGES": ADMIN_PAGES,
                "ALL_PAGES": ALL_PAGES,
            },
        ),
    )

    route_to_page(page)

    st.sidebar.markdown("---")
    st.sidebar.markdown("**Amprenta Multi-Omics Platform**")
    st.sidebar.markdown("Data stored in Postgres")
    st.sidebar.caption("Press `?` for keyboard shortcuts")
    if st.sidebar.button("🔄 Refresh Data"):
        st.cache_resource.clear()
        st.rerun()
    st.sidebar.markdown("---")
    st.sidebar.markdown("**API Access**")
    st.sidebar.markdown("[FastAPI Docs](http://localhost:8000/docs)")
    st.sidebar.markdown("[API Health](http://localhost:8000/health)")
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Quick Actions**")
    st.sidebar.markdown(
        """
- Ingest data via CLI scripts
- Use FastAPI for programmatic access
- Query Postgres directly for advanced queries
"""
    )


if __name__ == "__main__":
    main()

