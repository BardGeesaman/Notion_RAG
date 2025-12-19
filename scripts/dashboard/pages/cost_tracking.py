"""Cost Tracking page for managing project and experiment costs."""

from __future__ import annotations

from datetime import datetime, date, timedelta
from uuid import UUID

import pandas as pd
import plotly.express as px
import streamlit as st

from amprenta_rag.database.models import CostEntry, Project, Experiment
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_cost_tracking_page() -> None:
    """Render the Cost Tracking page."""
    st.header("💰 Cost Tracking")
    st.markdown("Track and analyze costs for projects and experiments.")

    tabs = st.tabs(["Overview", "Add Entry", "Entries"])

    with tabs[0]:
        _render_overview_tab()

    with tabs[1]:
        _render_add_entry_tab()

    with tabs[2]:
        _render_entries_tab()


def _render_overview_tab() -> None:
    """Render the overview tab with summary statistics and charts."""
    with db_session() as db:
        # Get all cost entries
        entries = db.query(CostEntry).all()

        if not entries:
            st.info("No cost entries yet. Add entries using the 'Add Entry' tab.")
            return

        # Calculate totals
        total_spend = sum(e.amount for e in entries)

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Entries", len(entries))
        with col2:
            st.metric("Total Spend", f"${total_spend:,.2f}")
        with col3:
            # Average per entry
            avg_spend = total_spend / len(entries) if entries else 0
            st.metric("Average per Entry", f"${avg_spend:,.2f}")

        st.markdown("---")

        # Spend by category (pie chart)
        category_data = {}
        for entry in entries:
            category_data[entry.category] = category_data.get(entry.category, 0) + entry.amount

        if category_data:
            st.subheader("Spend by Category")
            df_category = pd.DataFrame([
                {"Category": cat, "Amount": amt}
                for cat, amt in category_data.items()
            ])

            fig_pie = px.pie(
                df_category,
                values="Amount",
                names="Category",
                title="Cost Distribution by Category"
            )
            st.plotly_chart(fig_pie, use_container_width=True)

        # Spend by project (bar chart)
        project_data = {}
        for entry in entries:
            if entry.project_id:
                project = db.query(Project).filter(Project.id == entry.project_id).first()
                project_name = project.name if project else "Unknown"
            else:
                project_name = "Unassigned"
            project_data[project_name] = project_data.get(project_name, 0) + entry.amount

        if project_data:
            st.subheader("Spend by Project")
            df_project = pd.DataFrame([
                {"Project": proj, "Amount": amt}
                for proj, amt in sorted(project_data.items(), key=lambda x: x[1], reverse=True)
            ])

            fig_bar = px.bar(
                df_project,
                x="Project",
                y="Amount",
                title="Cost by Project",
                labels={"Amount": "Total Cost ($)", "Project": "Project Name"}
            )
            fig_bar.update_xaxes(tickangle=45)
            st.plotly_chart(fig_bar, use_container_width=True)

        # Recent entries summary
        st.markdown("---")
        st.subheader("Recent Entries")
        recent_entries = sorted(entries, key=lambda e: e.created_at, reverse=True)[:10]

        if recent_entries:
            recent_data = []
            for entry in recent_entries:
                project_name = ""
                if entry.project_id:
                    project = db.query(Project).filter(Project.id == entry.project_id).first()
                    project_name = project.name if project else "Unknown"

                recent_data.append({
                    "Date": entry.entry_date.strftime("%Y-%m-%d") if entry.entry_date else "",
                    "Category": entry.category,
                    "Description": entry.description[:50] + "..." if len(entry.description) > 50 else entry.description,
                    "Amount": f"${entry.amount:,.2f}",
                    "Project": project_name or "Unassigned",
                })

            df_recent = pd.DataFrame(recent_data)
            st.dataframe(df_recent, use_container_width=True, hide_index=True)


def _render_add_entry_tab() -> None:
    """Render the add entry form."""
    user = get_current_user()
    user_id = UUID(user.get("id")) if user and user.get("id") and user.get("id") != "test" else None

    if not user_id:
        st.warning("Please log in to add cost entries.")
        return

    with db_session() as db:
        st.subheader("Add New Cost Entry")

        with st.form("add_cost_entry", clear_on_submit=True):
            category = st.selectbox(
                "Category*",
                ["reagents", "equipment", "labor", "outsourcing", "other"]
            )

            amount = st.number_input("Amount*", min_value=0.0, step=0.01, format="%.2f")

            currency = st.selectbox("Currency", ["USD", "EUR", "GBP", "JPY"], index=0)

            description = st.text_area("Description*", placeholder="Describe this cost entry...", height=100)

            # Link to project or experiment
            link_type = st.radio("Link to", ["Project", "Experiment", "None"], horizontal=True)

            project_id = None
            experiment_id = None

            if link_type == "Project":
                projects = db.query(Project).order_by(Project.name).all()
                if projects:
                    project_options = {proj.name: proj.id for proj in projects}
                    selected_project_name = st.selectbox("Select Project", list(project_options.keys()))
                    project_id = project_options[selected_project_name]
                else:
                    st.info("No projects available.")

            elif link_type == "Experiment":
                experiments = db.query(Experiment).order_by(Experiment.name).limit(100).all()
                if experiments:
                    exp_options = {exp.name: exp.id for exp in experiments}
                    selected_exp_name = st.selectbox("Select Experiment", list(exp_options.keys()))
                    experiment_id = exp_options[selected_exp_name]
                else:
                    st.info("No experiments available.")

            entry_date = st.date_input("Entry Date*", value=date.today())

            submitted = st.form_submit_button("💾 Save Entry", type="primary")

            if submitted:
                if not description or not description.strip():
                    st.error("Description is required.")
                elif amount <= 0:
                    st.error("Amount must be greater than 0.")
                else:
                    try:
                        # Create entry
                        entry = CostEntry(
                            project_id=project_id,
                            experiment_id=experiment_id,
                            category=category,
                            description=description.strip(),
                            amount=amount,
                            currency=currency,
                            entry_date=datetime.combine(entry_date, datetime.min.time()),
                            created_by_id=user_id,
                        )

                        db.add(entry)
                        db.commit()
                        st.success(f"Cost entry added successfully!")
                        st.rerun()

                    except Exception as e:
                        st.error(f"Failed to add cost entry: {e}")


def _render_entries_tab() -> None:
    """Render the entries table with filters."""
    with db_session() as db:
        # Filters
        col1, col2, col3 = st.columns(3)

        with col1:
            # Category filter
            categories = db.query(CostEntry.category).distinct().all()
            category_options = ["All"] + [c[0] for c in categories if c[0]]
            selected_category = st.selectbox("Filter by Category", category_options)

        with col2:
            # Project filter
            projects_with_entries = (
                db.query(Project)
                .join(CostEntry, Project.id == CostEntry.project_id)
                .distinct()
                .all()
            )
            project_options = ["All"] + [proj.name for proj in projects_with_entries]
            selected_project = st.selectbox("Filter by Project", project_options)

        with col3:
            # Date range filter
            date_range = st.selectbox(
                "Date Range",
                ["All", "Last 7 days", "Last 30 days", "Last 90 days", "Last year"],
                index=0
            )

        # Build query
        query = db.query(CostEntry)

        if selected_category != "All":
            query = query.filter(CostEntry.category == selected_category)

        if selected_project != "All":
            proj = db.query(Project).filter(Project.name == selected_project).first()
            if proj:
                query = query.filter(CostEntry.project_id == proj.id)

        if date_range != "All":
            cutoff_date = datetime.now()
            if date_range == "Last 7 days":
                cutoff_date = cutoff_date - timedelta(days=7)
            elif date_range == "Last 30 days":
                cutoff_date = cutoff_date - timedelta(days=30)
            elif date_range == "Last 90 days":
                cutoff_date = cutoff_date - timedelta(days=90)
            elif date_range == "Last year":
                cutoff_date = cutoff_date - timedelta(days=365)

            query = query.filter(CostEntry.entry_date >= cutoff_date)

        entries = query.order_by(CostEntry.entry_date.desc()).all()

        st.metric("Filtered Entries", len(entries))

        if entries:
            # Calculate filtered total
            filtered_total = sum(e.amount for e in entries)
            st.metric("Filtered Total", f"${filtered_total:,.2f}")

            # Display entries table
            entry_data = []
            for entry in entries:
                project_name = ""
                experiment_name = ""

                if entry.project_id:
                    project = db.query(Project).filter(Project.id == entry.project_id).first()
                    project_name = project.name if project else "Unknown"

                if entry.experiment_id:
                    experiment = db.query(Experiment).filter(Experiment.id == entry.experiment_id).first()
                    experiment_name = experiment.name if experiment else "Unknown"

                entry_data.append({
                    "Date": entry.entry_date.strftime("%Y-%m-%d") if entry.entry_date else "",
                    "Category": entry.category,
                    "Description": entry.description,
                    "Amount": f"${entry.amount:,.2f} {entry.currency}",
                    "Project": project_name or "N/A",
                    "Experiment": experiment_name or "N/A",
                    "Created": entry.created_at.strftime("%Y-%m-%d %H:%M") if entry.created_at else "",
                })

            df_entries = pd.DataFrame(entry_data)
            st.dataframe(df_entries, use_container_width=True, hide_index=True)

            # Export button
            st.markdown("---")
            csv = df_entries.to_csv(index=False)
            st.download_button(
                label="📥 Download as CSV",
                data=csv,
                file_name=f"cost_entries_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv"
            )
        else:
            st.info("No entries found matching the selected filters.")


__all__ = ["render_cost_tracking_page"]

