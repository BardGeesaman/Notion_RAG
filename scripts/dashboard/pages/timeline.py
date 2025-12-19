"""Timeline page for viewing unified activity timeline."""
from __future__ import annotations

from datetime import datetime, timedelta

import streamlit as st

from amprenta_rag.utils.timeline import get_timeline
from scripts.dashboard.db_session import db_session


def render_timeline_page() -> None:
    """Render the Timeline page."""
    st.header("📅 Activity Timeline")
    st.markdown("View recent activity across experiments, compounds, samples, and discoveries.")

    with db_session() as db:
        # Filters
        col1, col2 = st.columns(2)
        with col1:
            entity_types = st.multiselect(
                "Filter by Type",
                ["experiment", "compound", "sample", "discovery"],
                default=["experiment", "compound", "sample", "discovery"],
                key="timeline_types"
            )
        with col2:
            days_back = st.selectbox(
                "Time Range",
                [1, 7, 30, 90, 365],
                index=2,
                format_func=lambda x: f"Last {x} days" if x < 365 else "All time",
                key="timeline_days"
            )

        # Get timeline
        timeline = get_timeline(limit=100, db=db)

        # Filter by selected types
        if entity_types:
            timeline = [item for item in timeline if item["type"] in entity_types]

        # Filter by date range
        if days_back < 365:
            cutoff = datetime.utcnow() - timedelta(days=days_back)
            timeline = [item for item in timeline if item["timestamp"] and item["timestamp"] >= cutoff]

        st.metric("Total Items", len(timeline))

        if timeline:
            st.markdown("---")
            st.subheader("Timeline")

            # Display timeline items
            for item in timeline:
                # Emoji icons per type
                icons = {
                    "experiment": "🧪",
                    "compound": "💊",
                    "sample": "🧫",
                    "discovery": "🔍",
                }
                icon = icons.get(item["type"], "📌")

                timestamp_str = item["timestamp"].strftime("%Y-%m-%d %H:%M") if item["timestamp"] else "Unknown"
                user_str = f" by {item['user']}" if item["user"] else ""

                col1, col2 = st.columns([4, 1])
                with col1:
                    st.markdown(f"{icon} **{item['name']}** ({item['type']}) - {timestamp_str}{user_str}")
                with col2:
                    # Determine which page to link to
                    page_map = {
                        "experiment": "Experiments",
                        "compound": "Chemistry",
                        "sample": "Sample Inventory",
                        "discovery": "Discovery Workflow",
                    }
                    page = page_map.get(item["type"], "Overview")
                    if st.button("View", key=f"view_{item['type']}_{item['id']}", use_container_width=True):
                        st.session_state["selected_page"] = page
                        st.session_state[f"selected_{item['type']}_id"] = item["id"]
                        st.rerun()

                st.divider()
        else:
            st.info("No activity found for the selected filters.")
