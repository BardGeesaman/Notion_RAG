"""Discovery Workflow - Automated repository scanning and import."""
import streamlit as st
from amprenta_rag.database.models import DiscoveryJob, DiscoveredStudy, HarvestSchedule
from amprenta_rag.database.session import db_session
from amprenta_rag.auth.session import get_current_user
from amprenta_rag.ingestion.harvest_scheduler import schedule_harvest


def render_discovery_workflow_page():
    st.title("🔍 Discovery Workflow")
    st.markdown("Scan repositories for new studies and import them")

    tab1, tab2, tab3, tab4 = st.tabs(["Run Discovery", "Pending Studies", "Job History", "Schedules"])

    with tab1:
        render_run_discovery_tab()

    with tab2:
        render_pending_tab()

    with tab3:
        render_history_tab()

    with tab4:
        render_schedules_tab()


def render_run_discovery_tab():
    st.subheader("Start New Discovery Job")
    user = get_current_user()

    col1, col2 = st.columns(2)
    with col1:
        repository = st.selectbox("Repository", ["GEO", "MetaboLights"])
    with col2:
        max_results = st.number_input("Max Results", min_value=10, max_value=200, value=50)

    query = st.text_input("Search Query", placeholder="e.g., ALS transcriptomics, Alzheimer proteomics")

    if st.button("🚀 Start Discovery", type="primary"):
        if not query:
            st.error("Please enter a search query")
            return

        with st.spinner(f"Scanning {repository}..."):
            from amprenta_rag.ingestion.discovery_service import run_discovery_job
            job_id = run_discovery_job(
                repository=repository,
                query=query,
                max_results=max_results,
                user_id=user.get("id") if user else None
            )
            st.success(f"Discovery job completed! Job ID: {job_id[:8]}...")
            st.rerun()


def render_pending_tab():
    st.subheader("Pending Studies")
    st.markdown("Review and import discovered studies")

    with db_session() as db:
        pending = db.query(DiscoveredStudy).filter(
            DiscoveredStudy.status == "new"
        ).order_by(DiscoveredStudy.discovered_at.desc()).limit(50).all()

        if not pending:
            st.info("No pending studies. Run a discovery job to find new studies!")
            return

        st.metric("Pending Review", len(pending))

        for study in pending:
            with st.expander(f"**{study.study_id}** ({study.repository}) - {study.omics_type or 'Unknown'}"):
                st.markdown(f"**Title:** {study.title or 'No title'}")
                st.markdown(f"**Organism:** {study.organism or 'Unknown'}")
                if study.description:
                    st.markdown(f"**Description:** {study.description[:300]}...")

                col1, col2, col3 = st.columns(3)
                with col1:
                    if st.button("✅ Import", key=f"import_{study.id}"):
                        from amprenta_rag.ingestion.discovery_service import import_discovered_study
                        exp_id = import_discovered_study(str(study.id))
                        if exp_id:
                            st.success(f"Imported as experiment!")
                            st.rerun()
                        else:
                            st.error("Import failed")
                with col2:
                    if st.button("⏭️ Skip", key=f"skip_{study.id}"):
                        study.status = "skipped"
                        db.commit()
                        st.rerun()
                with col3:
                    if st.button("🔗 View in Repo", key=f"view_{study.id}"):
                        if study.repository == "GEO":
                            st.markdown(f"[Open in GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={study.study_id})")
                        elif study.repository == "MetaboLights":
                            st.markdown(f"[Open in MetaboLights](https://www.ebi.ac.uk/metabolights/{study.study_id})")


def render_history_tab():
    st.subheader("Discovery Job History")

    with db_session() as db:
        jobs = db.query(DiscoveryJob).order_by(DiscoveryJob.created_at.desc()).limit(20).all()

        if not jobs:
            st.info("No discovery jobs yet")
            return

        for job in jobs:
            status_icon = {"completed": "✅", "failed": "❌", "running": "🔄", "pending": "⏳"}.get(job.status, "❓")
            with st.expander(f"{status_icon} {job.repository}: '{job.query}' ({job.status})"):
                col1, col2, col3 = st.columns(3)
                col1.metric("Found", job.studies_found)
                col2.metric("Imported", job.studies_imported)
                col3.metric("Status", job.status)

                if job.error_message:
                    st.error(f"Error: {job.error_message}")

                st.caption(f"Created: {job.created_at.strftime('%Y-%m-%d %H:%M') if job.created_at else 'Unknown'}")


def render_schedules_tab():
    st.subheader("Harvest Schedules")
    st.markdown("Manage automated repository harvesting schedules")

    user = get_current_user()

    with db_session() as db:
        # List active schedules
        st.markdown("### Active Schedules")
        schedules = db.query(HarvestSchedule).order_by(HarvestSchedule.created_at.desc()).all()

        if not schedules:
            st.info("No schedules yet. Create your first schedule below!")
        else:
            for schedule in schedules:
                status_icon = "🟢" if schedule.is_active else "🔴"
                with st.expander(f"{status_icon} **{schedule.name}** - {schedule.repository}"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown(f"**Query:** {schedule.query}")
                        st.markdown(f"**Interval:** Every {schedule.interval_hours} hours")
                    with col2:
                        st.markdown(f"**Last Run:** {schedule.last_run.strftime('%Y-%m-%d %H:%M') if schedule.last_run else 'Never'}")
                        st.markdown(f"**Next Run:** {schedule.next_run.strftime('%Y-%m-%d %H:%M') if schedule.next_run else 'Not scheduled'}")

                    col1, col2, col3 = st.columns(3)
                    with col1:
                        if st.button("🔄 Toggle Active", key=f"toggle_{schedule.id}"):
                            schedule.is_active = not schedule.is_active
                            db.commit()
                            if schedule.is_active:
                                schedule_harvest(str(schedule.id))
                            st.rerun()
                    with col2:
                        if st.button("🗑️ Delete", key=f"delete_{schedule.id}"):
                            db.delete(schedule)
                            db.commit()
                            st.success("Schedule deleted")
                            st.rerun()

        st.divider()

        # Create new schedule
        st.markdown("### Create New Schedule")
        with st.form("create_schedule"):
            name = st.text_input("Schedule Name*")
            repository = st.selectbox("Repository", ["GEO", "MetaboLights"])
            query = st.text_input("Search Query*", placeholder="e.g., ALS transcriptomics")
            interval_hours = st.number_input("Interval (hours)", min_value=1, max_value=168, value=24)

            submitted = st.form_submit_button("Create Schedule", type="primary")

            if submitted:
                if not name or not query:
                    st.error("Name and query are required")
                else:
                    new_schedule = HarvestSchedule(
                        name=name,
                        repository=repository,
                        query=query,
                        interval_hours=interval_hours,
                        is_active=True,
                        created_by_id=user.get("id") if user and user.get("id") != "test" else None,
                    )
                    db.add(new_schedule)
                    db.commit()

                    # Schedule the harvest job
                    schedule_harvest(str(new_schedule.id))

                    st.success(f"Schedule '{name}' created and activated!")
                    st.rerun()


if __name__ == "__main__":
    render_discovery_workflow_page()
