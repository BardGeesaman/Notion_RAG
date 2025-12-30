"""Job Queue Management dashboard."""

from __future__ import annotations

import os
from datetime import datetime
from typing import Any, Dict, Optional

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, params: Optional[Dict[str, Any]] = None, *, timeout: int = 60) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}", params=params or {})
    r.raise_for_status()
    return r.json()


def _api_post(path: str, json_body: Optional[Dict[str, Any]] = None, *, timeout: int = 60) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body or {})
    r.raise_for_status()
    return r.json()


def status_badge(status: str) -> str:
    """Return emoji badge for job status."""
    badges = {
        "pending": "🟡 Pending",
        "running": "🔵 Running", 
        "complete": "🟢 Complete",
        "completed": "🟢 Complete",
        "failed": "🔴 Failed",
        "cancelled": "⚪ Cancelled",
    }
    return badges.get(status, f"❓ {status}")


def format_duration(started_at: Optional[str], completed_at: Optional[str]) -> str:
    """Calculate and format job duration."""
    if not started_at:
        return "—"
    
    try:
        start = datetime.fromisoformat(started_at.replace('Z', '+00:00'))
        if completed_at:
            end = datetime.fromisoformat(completed_at.replace('Z', '+00:00'))
        else:
            end = datetime.now(start.tzinfo)
        
        duration = end - start
        total_seconds = int(duration.total_seconds())
        
        if total_seconds < 60:
            return f"{total_seconds}s"
        elif total_seconds < 3600:
            return f"{total_seconds // 60}m {total_seconds % 60}s"
        else:
            hours = total_seconds // 3600
            minutes = (total_seconds % 3600) // 60
            return f"{hours}h {minutes}m"
    except Exception:
        return "—"


def render_active_jobs_tab() -> None:
    """Render the Active Jobs tab."""
    st.subheader("🔄 Active Jobs")
    
    # Filter controls
    col1, col2, col3, col4 = st.columns([2, 2, 2, 1])
    
    with col1:
        job_type_filter = st.selectbox(
            "Job Type",
            options=["All", "genomics", "docking", "extraction", "sync", "single_cell"],
            index=0,
            key="active_job_type_filter"
        )
    
    with col2:
        status_filter = st.selectbox(
            "Status",
            options=["All", "pending", "running", "complete", "failed"],
            index=0,
            key="active_status_filter"
        )
    
    with col3:
        auto_refresh = st.checkbox("Auto-refresh (30s)", key="auto_refresh")
        if auto_refresh:
            st.rerun()
    
    with col4:
        if st.button("🔄 Refresh", key="refresh_active"):
            st.rerun()
    
    # Fetch jobs
    try:
        params = {"skip": 0, "limit": 50}
        if job_type_filter != "All":
            params["job_type"] = job_type_filter
        if status_filter != "All":
            params["status"] = status_filter
        
        response = _api_get("/api/v1/jobs", params=params)
        jobs = response.get("jobs", [])
        
        if not jobs:
            st.info("No jobs found matching the current filters.")
            return
        
        # Display jobs in a table
        job_data = []
        for job in jobs:
            # Calculate progress for certain job types
            progress = "—"
            if job["type"] == "genomics" and job.get("progress_percent") is not None:
                progress = f"{job['progress_percent']}%"
            elif job["type"] == "docking":
                total = job.get("total_compounds", 0)
                completed = job.get("completed_compounds", 0)
                if total > 0:
                    progress = f"{completed}/{total}"
            elif job["type"] == "extraction":
                total = job.get("file_count", 0)
                completed = job.get("completed_count", 0)
                if total > 0:
                    progress = f"{completed}/{total}"
            
            job_data.append({
                "Type": job["type"].title(),
                "Job ID": job["id"][:8] + "...",
                "Status": status_badge(job["status"]),
                "Progress": progress,
                "Started": job.get("started_at", "—"),
                "Duration": format_duration(job.get("started_at"), job.get("completed_at")),
                "Actions": job["id"],  # Will be used for action buttons
            })
        
        df = pd.DataFrame(job_data)
        
        # Display table
        st.dataframe(
            df.drop("Actions", axis=1),
            use_container_width=True,
            hide_index=True
        )
        
        # Action buttons for each job
        st.subheader("Job Actions")
        for i, job in enumerate(jobs):
            col1, col2, col3, col4 = st.columns([3, 2, 2, 2])
            
            with col1:
                st.text(f"{job['type'].title()} - {job['id'][:8]}...")
            
            with col2:
                if job["status"] in ["pending", "running"]:
                    if st.button("❌ Cancel", key=f"cancel_{job['id']}"):
                        try:
                            _api_post(f"/api/v1/jobs/{job['type']}/{job['id']}/cancel")
                            st.success("Job cancelled successfully!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to cancel job: {e}")
            
            with col3:
                if job["status"] == "failed":
                    if st.button("🔄 Retry", key=f"retry_{job['id']}"):
                        try:
                            _api_post(f"/api/v1/jobs/{job['type']}/{job['id']}/retry")
                            st.success("Job resubmitted successfully!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to retry job: {e}")
            
            with col4:
                if st.button("📋 Details", key=f"details_{job['id']}"):
                    try:
                        job_details = _api_get(f"/api/v1/jobs/{job['type']}/{job['id']}")
                        st.json(job_details)
                    except Exception as e:
                        st.error(f"Failed to fetch job details: {e}")
    
    except Exception as e:
        st.error(f"Failed to fetch jobs: {e}")


def render_job_history_tab() -> None:
    """Render the Job History tab."""
    st.subheader("📜 Job History")
    
    # Filter controls
    col1, col2, col3 = st.columns([2, 2, 2])
    
    with col1:
        st.selectbox(
            "Time Range",
            options=["Last 24 hours", "Last 7 days", "Last 30 days", "All time"],
            index=1,
            key="history_date_filter"
        )
    
    with col2:
        job_type_filter = st.selectbox(
            "Job Type",
            options=["All", "genomics", "docking", "extraction", "sync", "single_cell"],
            index=0,
            key="history_job_type_filter"
        )
    
    with col3:
        search_id = st.text_input("Search by Job ID", key="search_job_id")
    
    # Pagination controls
    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        page_size = st.selectbox("Page Size", options=[25, 50, 100], index=1, key="page_size")
    with col2:
        page_num = st.number_input("Page", min_value=1, value=1, key="page_num")
    
    # Fetch job history
    try:
        params = {
            "skip": (page_num - 1) * page_size,
            "limit": page_size
        }
        
        if job_type_filter != "All":
            params["job_type"] = job_type_filter
        
        response = _api_get("/api/v1/jobs", params=params)
        jobs = response.get("jobs", [])
        
        if search_id:
            jobs = [job for job in jobs if search_id.lower() in job["id"].lower()]
        
        if not jobs:
            st.info("No jobs found matching the current filters.")
            return
        
        # Display jobs with expandable details
        for job in jobs:
            with st.expander(f"{status_badge(job['status'])} {job['type'].title()} - {job['id'][:8]}... ({job.get('created_at', 'Unknown')})"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Job Details:**")
                    st.write(f"**ID:** {job['id']}")
                    st.write(f"**Type:** {job['type'].title()}")
                    st.write(f"**Status:** {status_badge(job['status'])}")
                    st.write(f"**Created:** {job.get('created_at', 'Unknown')}")
                    st.write(f"**Started:** {job.get('started_at', 'Not started')}")
                    st.write(f"**Completed:** {job.get('completed_at', 'Not completed')}")
                    st.write(f"**Duration:** {format_duration(job.get('started_at'), job.get('completed_at'))}")
                
                with col2:
                    st.write("**Type-specific Details:**")
                    if job["type"] == "genomics":
                        st.write(f"**Tool:** {job.get('tool', 'Unknown')}")
                        st.write(f"**Progress:** {job.get('progress_percent', 0)}%")
                        if job.get("error_message"):
                            st.error(f"**Error:** {job['error_message']}")
                    elif job["type"] == "docking":
                        st.write(f"**Compounds:** {job.get('completed_compounds', 0)}/{job.get('total_compounds', 0)}")
                        if job.get("error_log"):
                            st.error(f"**Error Log:** {job['error_log'][:200]}...")
                    elif job["type"] == "extraction":
                        st.write(f"**Files:** {job.get('completed_count', 0)}/{job.get('file_count', 0)}")
                    elif job["type"] == "sync":
                        st.write(f"**Source:** {job.get('source', 'Unknown')}")
                        st.write(f"**Records Synced:** {job.get('records_synced', 0)}")
                        st.write(f"**New Records:** {job.get('records_new', 0)}")
                    elif job["type"] == "single_cell":
                        st.write(f"**Cells:** {job.get('n_cells', 'Unknown')}")
                        st.write(f"**Genes:** {job.get('n_genes', 'Unknown')}")
        
        # Export button
        if st.button("📥 Export CSV", key="export_csv"):
            df = pd.DataFrame(jobs)
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name=f"job_history_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
    
    except Exception as e:
        st.error(f"Failed to fetch job history: {e}")


def render_statistics_tab() -> None:
    """Render the Statistics tab."""
    st.subheader("📊 Job Queue Statistics")
    
    try:
        # Fetch all jobs for statistics
        response = _api_get("/api/v1/jobs", params={"limit": 1000})
        jobs = response.get("jobs", [])
        
        if not jobs:
            st.info("No jobs available for statistics.")
            return
        
        # Calculate metrics
        total_jobs = len(jobs)
        running_jobs = len([j for j in jobs if j["status"] == "running"])
        failed_today = len([
            j for j in jobs 
            if j["status"] == "failed" and 
            j.get("created_at", "").startswith(datetime.now().strftime("%Y-%m-%d"))
        ])
        
        # Calculate average duration for completed jobs
        completed_jobs = [j for j in jobs if j["status"] in ["complete", "completed"] and j.get("started_at") and j.get("completed_at")]
        avg_duration = "—"
        if completed_jobs:
            total_duration = 0
            for job in completed_jobs:
                try:
                    start = datetime.fromisoformat(job["started_at"].replace('Z', '+00:00'))
                    end = datetime.fromisoformat(job["completed_at"].replace('Z', '+00:00'))
                    total_duration += (end - start).total_seconds()
                except Exception:
                    continue
            if total_duration > 0:
                avg_seconds = total_duration / len(completed_jobs)
                avg_duration = f"{int(avg_seconds // 60)}m {int(avg_seconds % 60)}s"
        
        # Metrics row
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Jobs", total_jobs)
        
        with col2:
            st.metric("Running", running_jobs)
        
        with col3:
            st.metric("Failed Today", failed_today)
        
        with col4:
            st.metric("Avg Duration", avg_duration)
        
        # Job type distribution
        st.subheader("Jobs by Type")
        type_counts = {}
        for job in jobs:
            job_type = job["type"]
            type_counts[job_type] = type_counts.get(job_type, 0) + 1
        
        if type_counts:
            df_types = pd.DataFrame(list(type_counts.items()), columns=["Type", "Count"])
            st.bar_chart(df_types.set_index("Type"))
        
        # Status distribution
        st.subheader("Jobs by Status")
        status_counts = {}
        for job in jobs:
            status = job["status"]
            status_counts[status] = status_counts.get(status, 0) + 1
        
        if status_counts:
            df_status = pd.DataFrame(list(status_counts.items()), columns=["Status", "Count"])
            st.bar_chart(df_status.set_index("Status"))
        
        # Queue information (placeholder - would need Celery inspection API)
        st.subheader("Queue Information")
        st.info("Queue depth information requires Celery inspection API integration.")
        
        # Flower dashboard link
        st.subheader("External Monitoring")
        st.markdown(
            "[🌸 Open Flower Dashboard](http://localhost:5555) - Real-time Celery monitoring",
            unsafe_allow_html=True
        )
        
    except Exception as e:
        st.error(f"Failed to fetch statistics: {e}")


def render_job_queue_page() -> None:
    """Render the Job Queue Management page."""
    from scripts.dashboard.auth import require_auth
    
    require_auth()
    
    st.set_page_config(page_title="Job Queue", page_icon="📋", layout="wide")
    st.title("📋 Job Queue Management")
    st.caption("Monitor and manage background job processing with Celery")
    
    # Three tabs
    tab1, tab2, tab3 = st.tabs(["🔄 Active Jobs", "📜 Job History", "📊 Statistics"])
    
    with tab1:
        render_active_jobs_tab()
    
    with tab2:
        render_job_history_tab()
    
    with tab3:
        render_statistics_tab()
    
    # Auto-refresh for active jobs tab
    if st.session_state.get("auto_refresh", False):
        import time
        time.sleep(30)
        st.rerun()
