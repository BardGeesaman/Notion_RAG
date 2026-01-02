"""Data Lifecycle Management Dashboard.

Provides UI for:
- Entity lifecycle status overview
- Quarantine queue management
- Bulk operations with preview
- Audit log viewing
"""

import os
from typing import Any, Dict, Optional

import httpx
import streamlit as st
from uuid import UUID

from scripts.dashboard.auth import require_auth


def _api_base() -> str:
    return os.environ.get("API_URL", "http://localhost:8000").rstrip("/")


def _api_get(path: str, params: Optional[Dict] = None) -> Any:
    try:
        with httpx.Client(timeout=60) as client:
            r = client.get(f"{_api_base()}{path}", params=params)
            r.raise_for_status()
            return r.json()
    except Exception:
        return None


def _api_post(path: str, payload: Dict) -> Any:
    try:
        with httpx.Client(timeout=60) as client:
            r = client.post(f"{_api_base()}{path}", json=payload)
            r.raise_for_status()
            return r.json()
    except Exception:
        return None


def render_data_lifecycle_page():
    """Main entry point for Data Lifecycle page."""
    require_auth()
    
    st.title("🗂️ Data Lifecycle Management")
    st.caption("Manage entity lifecycle status, quarantine, and bulk operations")
    
    # Tab layout
    tab_overview, tab_quarantine, tab_bulk, tab_audit = st.tabs([
        "📊 Overview",
        "🔒 Quarantine Queue",
        "⚡ Bulk Operations",
        "📜 Audit Log"
    ])
    
    with tab_overview:
        _render_overview_tab()
    
    with tab_quarantine:
        _render_quarantine_tab()
    
    with tab_bulk:
        _render_bulk_operations_tab()
    
    with tab_audit:
        _render_audit_tab()


def _render_overview_tab():
    """Overview of entity counts by lifecycle status."""
    st.subheader("Entity Status Distribution")
    
    # Fetch real stats from API
    stats = _api_get("/lifecycle/stats")
    
    if stats:
        # Display as metrics
        col1, col2, col3, col4 = st.columns(4)
        
        for col, (entity_type, counts) in zip(
            [col1, col2, col3, col4],
            [("Dataset", stats.get("dataset", {})),
             ("Experiment", stats.get("experiment", {})),
             ("Compound", stats.get("compound", {})),
             ("Signature", stats.get("signature", {}))]
        ):
            with col:
                st.markdown(f"**{entity_type}**")
                active = counts.get("active", 0)
                quarantined = counts.get("quarantined", 0)
                invalid = counts.get("invalid", 0)
                archived = counts.get("archived", 0)
                
                st.metric("Active", active)
                if quarantined > 0:
                    st.metric("Quarantined", quarantined, delta_color="off")
                if invalid > 0:
                    st.metric("Invalid", invalid, delta_color="inverse")
                if archived > 0:
                    st.metric("Archived", archived, delta_color="off")
    else:
        st.warning("Unable to fetch lifecycle stats")
    
    # Orphan Stats
    st.subheader("Orphaned Entities")
    orphan_stats = _api_get("/lifecycle/orphans")
    
    if orphan_stats:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Orphaned Features", orphan_stats.get("features", {}).get("count", 0))
        with col2:
            st.metric("Orphaned Signatures", orphan_stats.get("signatures", {}).get("count", 0))
        with col3:
            st.metric("Orphaned Embeddings", orphan_stats.get("embeddings", {}).get("count", 0))
        
        if st.button("🧹 Run Orphan Cleanup (Dry Run)", key="orphan_cleanup_preview"):
            # This would trigger the cleanup task - for now just show stats
            st.info("Orphan cleanup preview - no changes made")
    
    # Recent status changes
    st.subheader("Recent Status Changes")
    st.caption("Last 7 days of lifecycle transitions")
    
    # Placeholder for activity feed integration
    st.write("_Activity feed integration pending_")


def _render_quarantine_tab():
    """View and manage quarantined entities."""
    st.subheader("🔒 Quarantined Entities")
    
    entity_types = ["dataset", "experiment", "compound", "signature"]
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            entity_types,
            key="quarantine_entity_type"
        )
    
    with col2:
        st.info("Entities under review. Restore to make visible again, or archive to soft-delete.")
    
    # Action buttons for selected entities
    st.subheader("Quick Actions")
    
    col_id, col_action = st.columns([2, 1])
    
    with col_id:
        entity_id = st.text_input(
            "Entity ID (UUID)",
            placeholder="Enter UUID to check status or update",
            key="quarantine_entity_id"
        )
    
    with col_action:
        action = st.selectbox(
            "Action",
            ["Check Impact", "Restore (→ active)", "Archive (→ archived)", "Export Data"],
            key="quarantine_action"
        )
    
    if st.button("Execute", key="quarantine_execute", type="primary"):
        if not entity_id:
            st.error("Please enter an Entity ID")
        else:
            try:
                uuid_id = UUID(entity_id)
                
                if action == "Check Impact":
                    result = _api_get(f"/lifecycle/impact/{entity_type}/{uuid_id}")
                    if result:
                        st.json(result)
                    else:
                        st.error("Failed to fetch impact data")
                
                elif action == "Restore (→ active)":
                    result = _api_post("/lifecycle/status", {
                        "entity_type": entity_type,
                        "entity_id": str(uuid_id),
                        "new_status": "active",
                        "reason": "Restored from quarantine via dashboard"
                    })
                    if result and result.get("success"):
                        st.success(f"✅ {result.get('message')}")
                    else:
                        st.error(f"Failed: {result}")
                
                elif action == "Archive (→ archived)":
                    result = _api_post("/lifecycle/status", {
                        "entity_type": entity_type,
                        "entity_id": str(uuid_id),
                        "new_status": "archived",
                        "reason": "Archived from quarantine via dashboard"
                    })
                    if result and result.get("success"):
                        st.success(f"✅ {result.get('message')}")
                    else:
                        st.error(f"Failed: {result}")
                
                elif action == "Export Data":
                    result = _api_post("/lifecycle/export", {
                        "entity_type": entity_type,
                        "entity_id": str(uuid_id),
                    })
                    if result and not result.get("error"):
                        st.success(f"✅ Export ready - Checksum: {result.get('checksum', 'N/A')[:16]}...")
                        st.json(result.get("data", {}))
                        
                        # Download button
                        if st.button("📥 Download ZIP Package", key="download_export"):
                            st.info("Use POST /lifecycle/export/download endpoint directly for ZIP download")
                    else:
                        st.error(f"Export failed: {result.get('error', 'Unknown error')}")
                        
            except ValueError:
                st.error("Invalid UUID format")


def _render_bulk_operations_tab():
    """Bulk status updates and archive operations."""
    st.subheader("⚡ Bulk Operations")
    
    entity_types = ["dataset", "experiment", "compound", "signature"]
    status_options = ["active", "quarantined", "invalid", "archived"]
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            entity_types,
            key="bulk_entity_type"
        )
    
    with col2:
        new_status = st.selectbox(
            "Target Status",
            status_options,
            key="bulk_new_status"
        )
    
    # Entity IDs input
    entity_ids_input = st.text_area(
        "Entity IDs (one UUID per line, max 100)",
        height=150,
        placeholder="550e8400-e29b-41d4-a716-446655440000\n6ba7b810-9dad-11d1-80b4-00c04fd430c8\n...",
        key="bulk_entity_ids"
    )
    
    reason = st.text_input(
        "Reason for change",
        placeholder="Batch cleanup, data quality review, etc.",
        key="bulk_reason"
    )
    
    # Parse entity IDs
    entity_ids = []
    if entity_ids_input:
        lines = [line.strip() for line in entity_ids_input.split("\n") if line.strip()]
        for line in lines:
            try:
                entity_ids.append(str(UUID(line)))
            except ValueError:
                pass
    
    st.caption(f"Parsed {len(entity_ids)} valid UUIDs")
    
    # Preview button
    col_preview, col_execute = st.columns(2)
    
    with col_preview:
        if st.button("🔍 Preview Impact", key="bulk_preview"):
            if not entity_ids:
                st.warning("No valid UUIDs entered")
            else:
                result = _api_post("/lifecycle/bulk/preview", {
                    "entity_type": entity_type,
                    "entity_ids": entity_ids
                })
                if result:
                    st.session_state["bulk_preview_result"] = result
                    
                    # Display preview
                    st.subheader("Preview Results")
                    
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.metric("Entities", result.get("entity_count", 0))
                    with col_b:
                        can_proceed = result.get("can_proceed", False)
                        st.metric("Can Proceed", "✅ Yes" if can_proceed else "❌ No")
                    
                    if result.get("total_impact"):
                        st.write("**Impact Summary:**")
                        st.json(result["total_impact"])
                    
                    if result.get("blocking_entities"):
                        st.warning(f"⚠️ {len(result['blocking_entities'])} entities have blocking references")
                        with st.expander("View blocking details"):
                            st.json(result["blocking_entities"])
                else:
                    st.error("Failed to generate preview")
    
    with col_execute:
        # Only enable execute if preview was run
        preview_done = "bulk_preview_result" in st.session_state
        
        if st.button(
            "⚡ Execute Bulk Update",
            key="bulk_execute",
            type="primary",
            disabled=not preview_done
        ):
            if not entity_ids:
                st.warning("No valid UUIDs entered")
            elif not reason:
                st.warning("Please provide a reason for the change")
            else:
                # Use bulk/status for general updates, bulk/archive for archiving
                if new_status == "archived":
                    result = _api_post("/lifecycle/bulk/archive", {
                        "entity_type": entity_type,
                        "entity_ids": entity_ids,
                        "reason": reason,
                        "confirmed": True
                    })
                else:
                    result = _api_post("/lifecycle/bulk/status", {
                        "entity_type": entity_type,
                        "entity_ids": entity_ids,
                        "new_status": new_status,
                        "reason": reason
                    })
                
                if result:
                    success = result.get("success", 0)
                    failed = result.get("failed", 0)
                    total = result.get("total", 0)
                    
                    if failed == 0:
                        st.success(f"✅ All {success}/{total} entities updated successfully!")
                    else:
                        st.warning(f"⚠️ {success}/{total} succeeded, {failed} failed")
                        if result.get("errors"):
                            with st.expander("View errors"):
                                st.json(result["errors"])
                    
                    # Clear preview state
                    if "bulk_preview_result" in st.session_state:
                        del st.session_state["bulk_preview_result"]
                else:
                    st.error("Bulk operation failed")
    
    # After the existing bulk status update section, add:
    st.markdown("---")
    st.subheader("⚠️ Permanent Delete")
    st.warning("**DESTRUCTIVE ACTION** - Permanently deletes entities. Cannot be undone!")
    
    with st.expander("Permanent Delete (Advanced)", expanded=False):
        delete_entity_type = st.selectbox(
            "Entity Type",
            ["dataset", "experiment", "compound", "signature"],
            key="delete_entity_type"
        )
        
        delete_ids_input = st.text_area(
            "Entity IDs to DELETE (one per line)",
            height=100,
            placeholder="UUIDs of entities to permanently delete",
            key="delete_entity_ids"
        )
        
        delete_reason = st.text_input(
            "Reason for deletion",
            placeholder="GDPR request, data quality issue, etc.",
            key="delete_reason"
        )
        
        # Parse IDs
        delete_ids = []
        if delete_ids_input:
            for line in delete_ids_input.strip().split("\n"):
                try:
                    delete_ids.append(str(UUID(line.strip())))
                except ValueError:
                    pass
        
        st.caption(f"Parsed {len(delete_ids)} valid UUIDs")
        
        col_preview, col_delete = st.columns(2)
        
        with col_preview:
            if st.button("🔍 Preview Delete Impact", key="delete_preview"):
                if delete_ids:
                    result = _api_post("/lifecycle/bulk/preview", {
                        "entity_type": delete_entity_type,
                        "entity_ids": delete_ids
                    })
                    if result:
                        st.json(result)
        
        with col_delete:
            confirm_delete = st.checkbox("I understand this is permanent", key="confirm_permanent")
            if st.button("🗑️ PERMANENTLY DELETE", key="execute_delete", type="primary", disabled=not confirm_delete):
                if delete_ids and delete_reason:
                    result = _api_post("/lifecycle/bulk/delete", {
                        "entity_type": delete_entity_type,
                        "entity_ids": delete_ids,
                        "reason": delete_reason,
                        "confirmed": True,
                        "dry_run": False
                    })
                    if result:
                        st.success(f"Deleted {result.get('deleted', 0)} entities")
                        if result.get("errors"):
                            st.warning(f"Errors: {result['errors']}")
                else:
                    st.error("Provide IDs and reason")


def _render_audit_tab():
    """View lifecycle change audit log."""
    st.subheader("📜 Lifecycle Audit Log")
    
    # Filters
    col1, col2, col3 = st.columns(3)
    
    with col1:
        entity_type_filter = st.selectbox(
            "Entity Type",
            ["All", "dataset", "experiment", "compound", "signature"],
            key="audit_entity_type"
        )
    
    with col2:
        days_filter = st.selectbox(
            "Time Range",
            [7, 30, 90],
            format_func=lambda x: f"Last {x} days",
            key="audit_days"
        )
    
    with col3:
        st.write("")  # Spacer
        if st.button("🔄 Refresh", key="audit_refresh"):
            st.rerun()
    
    # Build query params
    params = f"?days={days_filter}&limit=50"
    if entity_type_filter != "All":
        params += f"&entity_type={entity_type_filter}"
    
    # Fetch audit data
    audit_data = _api_get(f"/lifecycle/audit{params}")
    
    if audit_data:
        total = audit_data.get("total", 0)
        entries = audit_data.get("entries", [])
        
        st.caption(f"Showing {len(entries)} of {total} entries")
        
        if entries:
            for entry in entries:
                with st.expander(
                    f"{entry.get('timestamp', 'N/A')[:19]} - {entry.get('entity_type', 'unknown')}:{entry.get('entity_id', 'N/A')[:8]}..."
                ):
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.write("**Old Value:**")
                        st.json(entry.get("old_value") or {})
                    with col_b:
                        st.write("**New Value:**")
                        st.json(entry.get("new_value") or {})
                    st.caption(f"Changed by: {entry.get('username', 'system')}")
        else:
            st.info("No lifecycle changes in this time range.")
    else:
        st.warning("Unable to fetch audit log")
