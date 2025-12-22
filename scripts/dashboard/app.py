"""
Streamlit dashboard entrypoint.

This is the canonical entry file for running the dashboard locally or in
containers (e.g., `streamlit run scripts/dashboard/app.py`).
"""

from __future__ import annotations

import streamlit as st

import scripts.dashboard.core as core


def main() -> None:
    """Render the Streamlit dashboard with sidebar navigation and dynamic routing."""
    st.set_page_config(page_title="Amprenta Dashboard", layout="wide")

    # Auth + session management
    user = core.check_authentication(core.AUTH_DISABLED)
    core.handle_session_timeout(core.AUTH_DISABLED, timeout_minutes=30)

    # Role-based page visibility (non-admin users don't see admin pages)
    visible_pages = list(core.ALL_PAGES)
    if not (user and user.get("role") == "admin"):
        visible_pages = [p for p in visible_pages if p not in set(core.ADMIN_PAGES)]

    page = core.render_sidebar(user, visible_pages, groups=core)
    core.route_to_page(page)


if __name__ == "__main__":
    main()


