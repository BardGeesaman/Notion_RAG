"""Schedule page for managing experiment scheduling and resource booking."""

from __future__ import annotations

from datetime import datetime, timedelta
from uuid import UUID

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from amprenta_rag.database.models import ScheduledEvent, Experiment
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_schedule_page() -> None:
    """Render the Schedule page."""
    st.header("📅 Experiment Schedule")
    st.markdown("Manage experiment scheduling and resource booking.")

    tabs = st.tabs(["Calendar", "Add Event"])

    with tabs[0]:
        _render_calendar_tab()

    with tabs[1]:
        _render_add_event_tab()


def _render_calendar_tab() -> None:
    """Render the calendar view with timeline/Gantt chart."""
    with db_session() as db:
        # Date range filter
        col1, col2 = st.columns(2)
        with col1:
            start_date = st.date_input("Start Date", value=datetime.now().date())
        with col2:
            end_date = st.date_input("End Date", value=(datetime.now() + timedelta(days=7)).date())

        # Resource filter
        resources = db.query(ScheduledEvent.resource_name).distinct().all()
        resource_options = ["All"] + [r[0] for r in resources if r[0]]
        selected_resource = st.selectbox("Filter by Resource", resource_options)

        # Build query
        query = db.query(ScheduledEvent).filter(
            ScheduledEvent.start_time >= datetime.combine(start_date, datetime.min.time()),
            ScheduledEvent.end_time <= datetime.combine(end_date, datetime.max.time())
        )

        if selected_resource != "All":
            query = query.filter(ScheduledEvent.resource_name == selected_resource)

        events = query.order_by(ScheduledEvent.start_time).all()

        st.metric("Scheduled Events", len(events))

        if events:
            # Prepare data for Gantt chart
            gantt_data = []
            for event in events:
                duration = (event.end_time - event.start_time).total_seconds() / 3600  # hours
                gantt_data.append({
                    "Resource": event.resource_name,
                    "Task": event.title,
                    "Start": event.start_time,
                    "End": event.end_time,
                    "Duration (hours)": duration,
                    "Event Type": event.event_type,
                    "Experiment": event.experiment.name if event.experiment else "N/A",
                })

            df = pd.DataFrame(gantt_data)

            # Create Gantt chart using Plotly
            fig = go.Figure()

            # Group by resource
            resources_list = df["Resource"].unique()
            colors = px.colors.qualitative.Set3

            for idx, resource in enumerate(resources_list):
                resource_df = df[df["Resource"] == resource]
                color = colors[idx % len(colors)]

                for _, row in resource_df.iterrows():
                    fig.add_trace(go.Scatter(
                        x=[row["Start"], row["End"]],
                        y=[resource],
                        mode='lines+markers',
                        name=row["Task"],
                        line=dict(width=20, color=color),
                        marker=dict(size=10),
                        text=f"{row['Task']}<br>{row['Start'].strftime('%Y-%m-%d %H:%M')} - {row['End'].strftime('%H:%M')}",
                        hovertemplate="<b>%{text}</b><br>Type: %{customdata[0]}<extra></extra>",
                        customdata=[[row["Event Type"]]],
                    ))

            fig.update_layout(
                title="Schedule Timeline",
                xaxis_title="Time",
                yaxis_title="Resource",
                height=400 + len(resources_list) * 50,
                showlegend=False,
                hovermode='closest',
            )

            st.plotly_chart(fig, use_container_width=True)

            # Events table
            st.markdown("---")
            st.subheader("Event Details")

            event_data = []
            for event in events:
                event_data.append({
                    "Title": event.title,
                    "Resource": event.resource_name,
                    "Type": event.event_type,
                    "Start": event.start_time.strftime("%Y-%m-%d %H:%M"),
                    "End": event.end_time.strftime("%Y-%m-%d %H:%M"),
                    "Experiment": event.experiment.name if event.experiment else "N/A",
                    "Created": event.created_at.strftime("%Y-%m-%d %H:%M") if event.created_at else "",
                })

            df_events = pd.DataFrame(event_data)
            st.dataframe(df_events, use_container_width=True, hide_index=True)

            # Detail view
            if events:
                selected_title = st.selectbox(
                    "Select event for details",
                    [e.title for e in events],
                    key="event_detail_select"
                )

                if selected_title:
                    event = next(e for e in events if e.title == selected_title)

                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**Title:** {event.title}")
                        st.write(f"**Event Type:** {event.event_type}")
                        st.write(f"**Resource:** {event.resource_name}")
                        st.write(f"**Start:** {event.start_time.strftime('%Y-%m-%d %H:%M')}")
                        st.write(f"**End:** {event.end_time.strftime('%Y-%m-%d %H:%M')}")

                    with col2:
                        if event.experiment:
                            st.write(f"**Experiment:** {event.experiment.name}")
                        if event.created_by_id:
                            from amprenta_rag.database.models import User
                            user = db.query(User).filter(User.id == event.created_by_id).first()
                            st.write(f"**Created By:** {user.username if user else 'Unknown'}")
                        st.write(f"**Created:** {event.created_at.strftime('%Y-%m-%d %H:%M') if event.created_at else 'N/A'}")

                    if event.notes:
                        st.markdown("---")
                        st.write("**Notes:**")
                        st.write(event.notes)
        else:
            st.info("No events scheduled for the selected date range.")


def _render_add_event_tab() -> None:
    """Render the add event form."""
    user = get_current_user()
    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None

    if not user_id:
        st.warning("Please log in to add scheduled events.")
        return

    with db_session() as db:
        st.subheader("Add New Scheduled Event")

        with st.form("add_scheduled_event", clear_on_submit=True):
            title = st.text_input("Title*", placeholder="e.g., Cell Culture Experiment - Day 1")

            event_type = st.selectbox(
                "Event Type*",
                ["experiment", "equipment", "meeting", "maintenance", "other"]
            )

            resource_name = st.text_input("Resource Name*", placeholder="e.g., Incubator A, Lab Room 101, Microscope B")

            col1, col2 = st.columns(2)
            with col1:
                start_time = st.datetime_input("Start Time*", value=datetime.now())
            with col2:
                end_time = st.datetime_input("End Time*", value=datetime.now() + timedelta(hours=2))

            # Link to experiment (optional)
            link_experiment = st.checkbox("Link to Experiment", value=False)
            experiment_id = None

            if link_experiment:
                experiments = db.query(Experiment).order_by(Experiment.name).all()
                if experiments:
                    exp_options = {exp.name: exp.id for exp in experiments}
                    selected_exp_name = st.selectbox("Select Experiment", list(exp_options.keys()))
                    experiment_id = exp_options[selected_exp_name]
                else:
                    st.info("No experiments available.")

            notes = st.text_area("Notes (optional)", placeholder="Additional notes about this event...", height=100)

            submitted = st.form_submit_button("💾 Save Event", type="primary")

            if submitted:
                if not title or not title.strip():
                    st.error("Title is required.")
                elif not resource_name or not resource_name.strip():
                    st.error("Resource name is required.")
                elif end_time <= start_time:
                    st.error("End time must be after start time.")
                else:
                    try:
                        # Check for conflicts
                        conflicts = db.query(ScheduledEvent).filter(
                            ScheduledEvent.resource_name == resource_name.strip(),
                            ScheduledEvent.start_time < end_time,
                            ScheduledEvent.end_time > start_time,
                        ).all()

                        if conflicts:
                            conflict_list = [f"{c.title} ({c.start_time.strftime('%Y-%m-%d %H:%M')} - {c.end_time.strftime('%H:%M')})" for c in conflicts]
                            st.warning("⚠️ Potential conflict detected with existing events:\n- " + "\n- ".join(conflict_list))

                            if not st.session_state.get("force_save_event", False):
                                if st.button("Force Save Anyway", key="force_save_btn"):
                                    st.session_state["force_save_event"] = True
                                    st.rerun()
                                return

                        # Create event
                        event = ScheduledEvent(
                            title=title.strip(),
                            event_type=event_type,
                            resource_name=resource_name.strip(),
                            start_time=start_time,
                            end_time=end_time,
                            experiment_id=experiment_id,
                            notes=notes.strip() if notes.strip() else None,
                            created_by_id=user_id,
                        )

                        db.add(event)
                        db.commit()
                        st.success(f"Event '{title}' scheduled successfully!")
                        st.session_state.pop("force_save_event", None)
                        st.rerun()

                    except Exception as e:
                        st.error(f"Failed to add event: {e}")


__all__ = ["render_schedule_page"]

