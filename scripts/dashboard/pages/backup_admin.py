"""Backup Administration page for managing database backups and project exports."""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Dict, List, Optional

import httpx
import pandas as pd
import streamlit as st

from scripts.dashboard.auth import require_admin


def render_backup_admin_page() -> None:
    """Render the Backup Administration page."""
    user = require_admin()
    if not user:
        return

    st.header("💾 Backup Administration")
    st.markdown("Manage database backups, scheduled tasks, and project exports.")

    tabs = st.tabs(["Manual Backup", "Scheduled Backups", "Project Export", "Restore Guide"])

    with tabs[0]:
        _render_manual_backup_tab()

    with tabs[1]:
        _render_scheduled_backups_tab()

    with tabs[2]:
        _render_project_export_tab()

    with tabs[3]:
        _render_restore_guide_tab()


def _render_manual_backup_tab() -> None:
    """Render the manual backup tab."""
    st.subheader("Manual Backup")
    st.markdown("Trigger immediate database backups and monitor progress.")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Create Backup")
        
        backup_type = st.selectbox(
            "Backup Type",
            options=["full", "incremental"],
            help="Full backup includes all data. Incremental includes only changes since last backup."
        )

        if st.button("🚀 Start Backup", type="primary", use_container_width=True):
            with st.spinner("Initiating backup..."):
                try:
                    result = _trigger_backup(backup_type)
                    if result:
                        st.success(f"✅ {result['message']}")
                        st.info(f"Task ID: `{result['task_id']}`")
                        
                        # Store task info in session state for tracking
                        st.session_state.backup_task_id = result['task_id']
                        st.session_state.backup_type = result['backup_type']
                        st.session_state.backup_started = datetime.now()
                        
                        st.rerun()
                    else:
                        st.error("❌ Failed to start backup")
                except Exception as e:
                    st.error(f"❌ Error: {str(e)}")

    with col2:
        st.markdown("### Current Status")
        
        # Show current task status if available
        if hasattr(st.session_state, 'backup_task_id'):
            st.info(f"**Task:** {st.session_state.backup_task_id[:8]}...")
            st.info(f"**Type:** {st.session_state.backup_type}")
            
            # Show elapsed time
            if hasattr(st.session_state, 'backup_started'):
                elapsed = datetime.now() - st.session_state.backup_started
                st.info(f"**Elapsed:** {str(elapsed).split('.')[0]}")
            
            if st.button("🔄 Check Status"):
                st.rerun()

    # Recent backups section
    st.markdown("### Recent Backups")
    _display_backup_history(limit=5)


def _render_scheduled_backups_tab() -> None:
    """Render the scheduled backups tab."""
    st.subheader("Scheduled Backups")
    st.markdown("Configure automatic backup schedules and view execution history.")

    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("### Schedule Configuration")
        
        # Mock schedule settings (in real implementation, these would be configurable)
        st.info("📅 **Daily Full Backup**\n\nScheduled: 1:00 AM UTC\nStatus: ✅ Enabled")
        st.info("🧹 **Weekly Cleanup**\n\nScheduled: Sunday 4:00 AM UTC\nStatus: ✅ Enabled")
        st.info("🔍 **Weekly Verification**\n\nScheduled: Sunday 5:00 AM UTC\nStatus: ✅ Enabled")
        
        # Toggle switches (mock functionality)
        st.checkbox("Enable Daily Backups", value=True)
        st.checkbox("Enable Weekly Cleanup", value=True)
        st.checkbox("Enable Weekly Verification", value=True)
        
        if st.button("💾 Save Schedule Settings"):
            st.success("✅ Schedule settings saved")

    with col2:
        st.markdown("### Next Scheduled Runs")
        
        # Calculate next run times (mock)
        now = datetime.now()
        next_daily = now.replace(hour=1, minute=0, second=0, microsecond=0)
        if next_daily <= now:
            next_daily += timedelta(days=1)
        
        next_cleanup = now + timedelta(days=(6 - now.weekday()))
        next_cleanup = next_cleanup.replace(hour=4, minute=0, second=0, microsecond=0)
        
        next_verify = next_cleanup.replace(hour=5)
        
        st.metric("Next Daily Backup", next_daily.strftime("%Y-%m-%d %H:%M UTC"))
        st.metric("Next Cleanup", next_cleanup.strftime("%Y-%m-%d %H:%M UTC"))
        st.metric("Next Verification", next_verify.strftime("%Y-%m-%d %H:%M UTC"))

    # Scheduled backup history
    st.markdown("### Execution History")
    _display_backup_history(limit=10)


def _render_project_export_tab() -> None:
    """Render the project export tab."""
    st.subheader("Project Export")
    st.markdown("Export selected project data as downloadable ZIP packages.")

    # Entity selection
    st.markdown("### Select Entities to Export")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**Programs**")
        program_ids = st.text_area(
            "Program UUIDs (one per line)",
            height=100,
            help="Enter program UUIDs to export, one per line"
        )
    
    with col2:
        st.markdown("**Experiments**")
        experiment_ids = st.text_area(
            "Experiment UUIDs (one per line)",
            height=100,
            help="Enter experiment UUIDs to export, one per line"
        )
    
    with col3:
        st.markdown("**Compounds**")
        compound_ids = st.text_area(
            "Compound UUIDs (one per line)",
            height=100,
            help="Enter compound UUIDs to export, one per line"
        )

    # Export options
    include_related = st.checkbox(
        "Include Related Entities",
        value=True,
        help="Automatically include datasets, features, and signatures related to selected entities"
    )

    # Export button
    if st.button("📦 Generate Export", type="primary", use_container_width=True):
        # Parse UUIDs
        programs = [uid.strip() for uid in program_ids.split('\n') if uid.strip()]
        experiments = [uid.strip() for uid in experiment_ids.split('\n') if uid.strip()]
        compounds = [uid.strip() for uid in compound_ids.split('\n') if uid.strip()]
        
        if not any([programs, experiments, compounds]):
            st.error("❌ Please specify at least one entity to export")
            return
        
        with st.spinner("Generating export package..."):
            try:
                result = _create_export(programs, experiments, compounds, include_related)
                if result:
                    st.success(f"✅ {result['message']}")
                    st.info(f"Export size: {result['export_size_bytes']} bytes")
                    st.info(f"Entities: {result['entities']}")
                    
                    # In a real implementation, provide download link
                    st.markdown("**Download will be available once export storage is implemented**")
                else:
                    st.error("❌ Failed to create export")
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")

    # Export history
    st.markdown("### Export History")
    st.info("Export history will be displayed here once export storage is implemented")


def _render_restore_guide_tab() -> None:
    """Render the restore guide tab."""
    st.subheader("Disaster Recovery Guide")
    st.markdown("Step-by-step instructions for restoring from backups.")

    # Emergency contacts
    with st.expander("🚨 Emergency Contacts", expanded=True):
        st.markdown("""
        **Database Administrator:** admin@company.com
        **DevOps Team:** devops@company.com
        **Emergency Hotline:** +1-555-BACKUP
        """)

    # Restore procedures
    st.markdown("### Restore Procedures")
    
    tab1, tab2, tab3 = st.tabs(["Full Restore", "Partial Restore", "Point-in-Time"])
    
    with tab1:
        st.markdown("""
        #### Full Database Restore
        
        **Prerequisites:**
        - Database server access
        - Latest backup file
        - PostgreSQL admin credentials
        
        **Steps:**
        1. **Stop application services**
           ```bash
           docker-compose down
           ```
        
        2. **Download backup file**
           - Use the backup download endpoint
           - Verify checksum integrity
        
        3. **Restore database**
           ```bash
           # Drop existing database (CAUTION!)
           dropdb amprenta_production
           
           # Create new database
           createdb amprenta_production
           
           # Restore from backup
           gunzip -c backup_full_YYYYMMDD_HHMMSS.sql.gz | psql amprenta_production
           ```
        
        4. **Verify restore**
           ```bash
           psql amprenta_production -c "SELECT COUNT(*) FROM programs;"
           ```
        
        5. **Restart services**
           ```bash
           docker-compose up -d
           ```
        """)
    
    with tab2:
        st.markdown("""
        #### Partial Data Restore
        
        **Use Cases:**
        - Recover specific tables
        - Restore deleted records
        - Fix data corruption
        
        **Steps:**
        1. **Extract specific tables**
           ```bash
           pg_restore -t programs -t experiments backup.sql.gz
           ```
        
        2. **Selective import**
           ```bash
           psql amprenta_production -f extracted_data.sql
           ```
        
        3. **Verify data integrity**
           - Check foreign key constraints
           - Validate data relationships
        """)
    
    with tab3:
        st.markdown("""
        #### Point-in-Time Recovery
        
        **Requirements:**
        - WAL archiving enabled
        - Base backup + WAL files
        
        **Process:**
        1. **Restore base backup**
        2. **Apply WAL files up to target time**
        3. **Set recovery target**
           ```
           recovery_target_time = '2024-01-15 14:30:00'
           ```
        """)

    # Verification checklist
    st.markdown("### Post-Restore Verification")
    
    checklist_items = [
        "Database connectivity established",
        "All tables present and populated",
        "Foreign key constraints valid",
        "Application services running",
        "User authentication working",
        "API endpoints responding",
        "Data integrity checks passed",
        "Performance benchmarks met"
    ]
    
    st.markdown("**Verification Checklist:**")
    for item in checklist_items:
        st.checkbox(item, key=f"checklist_{item.replace(' ', '_')}")


def _trigger_backup(backup_type: str) -> Optional[Dict]:
    """Trigger a manual backup via API."""
    try:
        # Use mock user ID for testing when auth is disabled
        headers = {"X-User-Id": "00000000-0000-0000-0000-000000000001"}
        
        with httpx.Client() as client:
            response = client.post(
                "http://localhost:8000/api/v1/backup/database",
                json={"backup_type": backup_type},
                headers=headers,
                timeout=30.0
            )
            
            if response.status_code == 202:
                return response.json()
            else:
                st.error(f"API Error: {response.status_code} - {response.text}")
                return None
                
    except Exception as e:
        st.error(f"Connection error: {str(e)}")
        return None


def _create_export(programs: List[str], experiments: List[str], compounds: List[str], include_related: bool) -> Optional[Dict]:
    """Create a project export via API."""
    try:
        # Use mock user ID for testing when auth is disabled
        headers = {"X-User-Id": "00000000-0000-0000-0000-000000000001"}
        
        # Prepare request data
        export_data = {
            "include_related": include_related
        }
        
        if programs:
            export_data["program_ids"] = programs
        if experiments:
            export_data["experiment_ids"] = experiments
        if compounds:
            export_data["compound_ids"] = compounds
        
        with httpx.Client() as client:
            response = client.post(
                "http://localhost:8000/api/v1/backup/export",
                json=export_data,
                headers=headers,
                timeout=60.0
            )
            
            if response.status_code == 202:
                return response.json()
            else:
                st.error(f"API Error: {response.status_code} - {response.text}")
                return None
                
    except Exception as e:
        st.error(f"Connection error: {str(e)}")
        return None


def _display_backup_history(limit: int = 10) -> None:
    """Display backup history table."""
    try:
        # Use mock user ID for testing when auth is disabled
        headers = {"X-User-Id": "00000000-0000-0000-0000-000000000001"}
        
        with httpx.Client() as client:
            response = client.get(
                f"http://localhost:8000/api/v1/backup/history?per_page={limit}",
                headers=headers,
                timeout=30.0
            )
            
            if response.status_code == 200:
                data = response.json()
                backups = data.get("items", [])
                
                if backups:
                    # Convert to DataFrame for display
                    df_data = []
                    for backup in backups:
                        df_data.append({
                            "ID": backup["id"][:8] + "...",
                            "Type": backup["backup_type"],
                            "Status": backup["status"],
                            "Size": f"{backup.get('file_size_bytes', 0) // 1024 // 1024} MB" if backup.get('file_size_bytes') else "N/A",
                            "Created": backup["created_at"][:19].replace("T", " "),
                            "Completed": backup.get("completed_at", "")[:19].replace("T", " ") if backup.get("completed_at") else "N/A"
                        })
                    
                    df = pd.DataFrame(df_data)
                    st.dataframe(df, use_container_width=True)
                    
                    # Show total count
                    st.caption(f"Showing {len(backups)} of {data.get('total', 0)} total backups")
                else:
                    st.info("No backup history available")
            else:
                st.error(f"Failed to load backup history: {response.status_code}")
                
    except Exception as e:
        st.error(f"Error loading backup history: {str(e)}")
