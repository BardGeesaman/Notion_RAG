"""
Activity Feed page for the dashboard.

Shows a timeline view of all activity events with filtering options.
"""

from __future__ import annotations

import os
from datetime import datetime, timedelta
from typing import List, Optional
import httpx
import streamlit as st

API_BASE = os.environ.get("API_URL", "http://localhost:8000")
ACTIVITY_ENDPOINT = f"{API_BASE}/api/v1/activity/feed"
PROGRAMS_ENDPOINT = f"{API_BASE}/api/v1/programs"


@st.cache_data(ttl=60)
def fetch_programs():
    """Fetch available programs for filtering."""
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(PROGRAMS_ENDPOINT)
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load programs: {e}")
        return []


@st.cache_data(ttl=30)
def fetch_activity_feed(
    program_id: Optional[str] = None,
    event_types: Optional[List[str]] = None,
    since: Optional[str] = None,
    limit: int = 50,
    offset: int = 0,
):
    """Fetch activity feed with filters."""
    params = {
        "limit": limit,
        "offset": offset,
    }
    
    if program_id:
        params["program_id"] = program_id
    
    if event_types:
        # For multiple event types, we'll filter client-side for simplicity
        # In a real implementation, the API should support multiple event_type filters
        pass
    
    if since:
        params["since"] = since
    
    try:
        with httpx.Client(timeout=10) as client:
            resp = client.get(ACTIVITY_ENDPOINT, params=params)
            resp.raise_for_status()
            events = resp.json()
            
            # Client-side filtering for event types if specified
            if event_types:
                events = [e for e in events if e.get("event_type") in event_types]
            
            return events
    except httpx.HTTPError as e:  # noqa: BLE001
        st.error(f"Failed to load activity feed: {e}")
        return []


def get_event_icon(event_type: str) -> str:
    """Get icon for activity event type."""
    icon_map = {
        "compound_added": "🧪",
        "experiment_created": "📊",
        "model_trained": "🤖",
        "hit_confirmed": "🎯",
        "status_changed": "📝",
        "notebook_reviewed": "📖",
    }
    return icon_map.get(event_type, "📋")


def format_event_description(event: dict) -> str:
    """Format event description for display."""
    event_type = event.get("event_type", "unknown")
    target_name = event.get("target_name", "Unknown")
    target_type = event.get("target_type", "item")
    
    action_map = {
        "compound_added": f"added compound **{target_name}**",
        "experiment_created": f"created experiment **{target_name}**",
        "model_trained": f"trained model **{target_name}**",
        "hit_confirmed": f"confirmed hit **{target_name}**",
        "status_changed": f"updated status of **{target_name}**",
        "notebook_reviewed": f"reviewed notebook **{target_name}**",
    }
    
    return action_map.get(event_type, f"performed {event_type} on **{target_name}**")


def format_timestamp(created_at: str) -> str:
    """Format timestamp for display."""
    try:
        dt = datetime.fromisoformat(created_at.replace('Z', '+00:00'))
        return dt.strftime("%Y-%m-%d %H:%M")
    except Exception:
        return created_at


def render_activity_card(event: dict):
    """Render a single activity event card."""
    event_type = event.get("event_type", "unknown")
    actor_id = event.get("actor_id")
    target_id = event.get("target_id")
    created_at = event.get("created_at", "")
    metadata = event.get("metadata", {})
    
    # Card container
    with st.container():
        col1, col2, col3 = st.columns([0.5, 4, 1.5])
        
        with col1:
            # Event icon
            icon = get_event_icon(event_type)
            st.markdown(f"**{icon}**")
        
        with col2:
            # Actor and action
            actor_name = metadata.get("actor_name", "System")
            if actor_id:
                st.markdown(f"**{actor_name}**")
            else:
                st.markdown("**System**")
            
            # Action description
            description = format_event_description(event)
            st.markdown(description)
            
            # Additional metadata if available
            if metadata:
                details = []
                if "program_name" in metadata:
                    details.append(f"Program: {metadata['program_name']}")
                if "experiment_type" in metadata:
                    details.append(f"Type: {metadata['experiment_type']}")
                if details:
                    st.caption(" • ".join(details))
        
        with col3:
            # Timestamp
            timestamp = format_timestamp(created_at)
            st.caption(timestamp)
            
            # Target link (if applicable)
            if target_id and event.get("target_type") in ["compound", "experiment", "dataset"]:
                if st.button("View", key=f"view_{event['id']}", help="View target"):
                    # Navigate to appropriate page based on target_type
                    target_type = event.get("target_type")
                    if target_type == "compound":
                        st.query_params["compound_id"] = str(target_id)
                        st.switch_page("pages/compound_details.py")
                    elif target_type == "experiment":
                        st.query_params["experiment_id"] = str(target_id)
                        st.switch_page("pages/experiment_details.py")
                    elif target_type == "dataset":
                        st.query_params["dataset_id"] = str(target_id)
                        st.switch_page("pages/dataset_details.py")
        
        st.markdown("---")


def render_activity_feed_page():
    """Render the main activity feed page."""
    st.title("📋 Activity Feed")
    st.markdown("Track all activities across your research programs.")
    
    # Sidebar filters
    with st.sidebar:
        st.header("Filters")
        
        # Program filter
        programs = fetch_programs()
        program_options = {"All Programs": None}
        for program in programs:
            program_options[program.get("name", "Unknown")] = program.get("id")
        
        selected_program_name = st.selectbox(
            "Program",
            options=list(program_options.keys()),
            index=0
        )
        selected_program_id = program_options[selected_program_name]
        
        # Event type filter
        available_event_types = [
            "compound_added",
            "experiment_created",
            "model_trained",
            "hit_confirmed",
            "status_changed",
            "notebook_reviewed",
        ]
        
        selected_event_types = st.multiselect(
            "Event Types",
            options=available_event_types,
            default=[],
            format_func=lambda x: x.replace("_", " ").title()
        )
        
        # Date range filter
        date_options = {
            "All Time": None,
            "Last 24 Hours": (datetime.now() - timedelta(days=1)).isoformat(),
            "Last Week": (datetime.now() - timedelta(weeks=1)).isoformat(),
            "Last Month": (datetime.now() - timedelta(days=30)).isoformat(),
        }
        
        selected_date_name = st.selectbox(
            "Time Range",
            options=list(date_options.keys()),
            index=1  # Default to Last 24 Hours
        )
        since_date = date_options[selected_date_name]
        
        # Custom date range
        if selected_date_name == "All Time":
            custom_since = st.date_input(
                "Custom Start Date (Optional)",
                value=None
            )
            if custom_since:
                since_date = custom_since.isoformat()
    
    # Main content
    col1, col2 = st.columns([3, 1])
    
    with col1:
        st.markdown("### Recent Activity")
    
    with col2:
        # Pagination controls
        if "activity_page" not in st.session_state:
            st.session_state.activity_page = 0
        
        page = st.session_state.activity_page
        limit = 20
        offset = page * limit
    
    # Fetch and display activity feed
    events = fetch_activity_feed(
        program_id=selected_program_id,
        event_types=selected_event_types if selected_event_types else None,
        since=since_date,
        limit=limit,
        offset=offset,
    )
    
    if not events:
        st.info("No activity found matching the selected filters.")
    else:
        # Display events
        for event in events:
            render_activity_card(event)
        
        # Pagination controls
        col1, col2, col3 = st.columns([1, 2, 1])
        
        with col1:
            if page > 0:
                if st.button("← Previous", use_container_width=True):
                    st.session_state.activity_page = page - 1
                    st.rerun()
        
        with col2:
            st.markdown(f"<center>Page {page + 1}</center>", unsafe_allow_html=True)
        
        with col3:
            if len(events) == limit:  # Might have more pages
                if st.button("Next →", use_container_width=True):
                    st.session_state.activity_page = page + 1
                    st.rerun()
    
    # Summary statistics
    st.markdown("---")
    with st.expander("📊 Activity Summary", expanded=False):
        if events:
            # Count events by type
            event_type_counts = {}
            for event in events:
                event_type = event.get("event_type", "unknown")
                event_type_counts[event_type] = event_type_counts.get(event_type, 0) + 1
            
            st.markdown("**Event Types in Current View:**")
            for event_type, count in sorted(event_type_counts.items()):
                icon = get_event_icon(event_type)
                st.markdown(f"- {icon} {event_type.replace('_', ' ').title()}: {count}")
        else:
            st.markdown("No events to summarize.")


if __name__ == "__main__":
    render_activity_feed_page()
