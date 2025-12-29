"""Scientist's Cockpit - Overview Dashboard."""

from __future__ import annotations

from datetime import datetime

import streamlit as st

from scripts.dashboard.components.cockpit_widgets import (
    render_activity_widget,
    render_alerts_widget,
    render_shortcuts_widget,
    render_stats_widget,
    render_tasks_widget,
)


def render_cockpit_page() -> None:
    """Render the Scientist's Cockpit overview dashboard."""
    from scripts.dashboard.auth import require_auth
    user = require_auth()
    
    # Welcome banner
    st.title("🚀 Scientist's Cockpit")
    username = user.get("username", "Scientist") if user else "Scientist"
    st.caption(f"Welcome, {username} • {datetime.now().strftime('%B %d, %Y')}")
    
    # Stats row (single widget with internal columns)
    render_stats_widget()
    
    st.markdown("---")
    
    # Activity and Alerts row
    col1, col2 = st.columns(2)
    
    with col1:
        render_activity_widget()
    
    with col2:
        render_alerts_widget()
    
    st.markdown("---")
    
    # Tasks and Shortcuts row
    col1, col2 = st.columns(2)
    
    with col1:
        render_tasks_widget()
    
    with col2:
        render_shortcuts_widget()


if __name__ == "__main__":
    render_cockpit_page()

