"""Discovery Workflow - Automated repository scanning and import."""
import streamlit as st
from datetime import datetime
from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import DiscoveryJob, DiscoveredStudy
from amprenta_rag.auth.session import get_current_user


def render_discovery_workflow_page():
    st.title("🔍 Discovery Workflow")
    st.markdown("Scan repositories for new studies and import them")

    tab1, tab2, tab3 = st.tabs(["Run Discovery", "Pending Studies", "Job History"])

    with tab1:
        render_run_discovery_tab()

    with tab2:
        render_pending_tab()

    with tab3:
        render_history_tab()


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

    db_gen = get_db()
    db = next(db_gen)
    try:
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
    finally:
        db_gen.close()


def render_history_tab():
    st.subheader("Discovery Job History")

    db_gen = get_db()
    db = next(db_gen)
    try:
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
    finally:
        db_gen.close()


if __name__ == "__main__":
    render_discovery_workflow_page()
