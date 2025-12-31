"""Provenance Ledger page."""

from __future__ import annotations

import os
from typing import Any
from uuid import UUID

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: dict, *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_provenance_ledger_page() -> None:
    """Render the Provenance Ledger page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📜 Provenance Ledger")
    st.caption("Browse entity version history, compare versions, and restore previous states.")
    
    tab1, tab2, tab3 = st.tabs(["Version History", "Compare Versions", "Restore Version"])
    
    with tab1:
        render_version_history_tab()
    
    with tab2:
        render_compare_versions_tab()
    
    with tab3:
        render_restore_version_tab()


def render_version_history_tab() -> None:
    """Render version history tab."""
    st.subheader("Version History")
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment"],
            key="history_entity_type",
        )
    
    with col2:
        entity_id = st.text_input(
            "Entity ID", 
            placeholder="UUID",
            key="history_entity_id",
        )
    
    if st.button("Load Versions", type="primary", disabled=not entity_id):
        try:
            # Validate UUID format
            UUID(entity_id)
            
            versions = _api_get(f"/api/v1/versions/{entity_type.lower()}/{entity_id}")
            st.session_state["version_history"] = versions
            st.session_state["history_entity_type"] = entity_type
            st.session_state["history_entity_id"] = entity_id
            
        except ValueError:
            st.error("Invalid UUID format")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                st.error("Entity not found or has no versions")
            else:
                st.error(f"Failed to load versions: {e}")
        except Exception as e:
            st.error(f"Failed to load versions: {e}")
    
    # Display versions
    versions = st.session_state.get("version_history")
    
    if versions:
        if versions:
            # Create dataframe with key columns
            df_data = []
            for version in versions:
                df_data.append({
                    "ID": version["id"],
                    "Version": version["version_number"],
                    "Checksum": version["checksum_sha256"][:16] + "...",
                    "Created": version["created_at"],
                    "Summary": version["change_summary"],
                })
            
            df = pd.DataFrame(df_data)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.caption(f"Total versions: {len(versions)}")
            
            # Expandable version details
            st.markdown("### Version Details")
            for version in versions:
                with st.expander(f"Version {version['version_number']} - {version['change_summary']}"):
                    st.json(version["data_snapshot"])
        else:
            st.info("No versions found for this entity")


def render_compare_versions_tab() -> None:
    """Render compare versions tab."""
    st.subheader("Compare Versions")
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment"],
            key="compare_entity_type",
        )
    
    with col2:
        entity_id = st.text_input(
            "Entity ID", 
            placeholder="UUID",
            key="compare_entity_id",
        )
    
    col3, col4 = st.columns(2)
    
    with col3:
        version_a = st.number_input(
            "Version A", 
            min_value=1, 
            value=1,
            key="version_a",
        )
    
    with col4:
        version_b = st.number_input(
            "Version B", 
            min_value=1, 
            value=2,
            key="version_b",
        )
    
    if st.button("Compare", type="primary", disabled=not entity_id):
        try:
            # Validate UUID format
            UUID(entity_id)
            
            comparison = _api_post(
                f"/api/v1/versions/{entity_type.lower()}/{entity_id}/compare",
                {
                    "version_a": version_a,
                    "version_b": version_b,
                }
            )
            st.session_state["version_comparison"] = comparison
            
        except ValueError:
            st.error("Invalid UUID format")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                st.error("Entity or versions not found")
            else:
                st.error(f"Failed to compare versions: {e}")
        except Exception as e:
            st.error(f"Failed to compare versions: {e}")
    
    # Display comparison
    comparison = st.session_state.get("version_comparison")
    
    if comparison:
        diff = comparison.get("diff", {})
        
        # Added fields
        added = diff.get("added", {})
        if added:
            st.markdown("### ✅ Added Fields")
            for key, value in added.items():
                st.success(f"**{key}**: `{value}`")
        
        # Removed fields
        removed = diff.get("removed", {})
        if removed:
            st.markdown("### ❌ Removed Fields")
            for key, value in removed.items():
                st.error(f"**{key}**: `{value}`")
        
        # Changed fields
        changed = diff.get("changed", {})
        if changed:
            st.markdown("### 🔄 Changed Fields")
            for key, changes in changed.items():
                st.warning(f"**{key}**: `{changes['old']}` → `{changes['new']}`")
        
        if not added and not removed and not changed:
            st.info("No differences found between versions")


def render_restore_version_tab() -> None:
    """Render restore version tab."""
    st.subheader("Restore Version")
    
    # Admin check
    user_role = st.session_state.get("user", {}).get("role")
    if user_role != "admin":
        st.error("🔒 Admin access required for version restore")
        return
    
    st.warning("⚠️ **Warning**: Restoring a version will create a new version with the old data. This action is logged in the audit trail.")
    
    version_id = st.text_input(
        "Version ID",
        placeholder="UUID of version to restore",
        key="restore_version_id",
    )
    
    reason = st.text_input(
        "Reason",
        placeholder="Why are you restoring this version?",
        key="restore_reason",
    )
    
    confirm = st.checkbox(
        "I confirm this restore action",
        key="restore_confirm",
    )
    
    restore_button = st.button(
        "Restore", 
        type="primary",
        disabled=not (version_id and reason and confirm),
    )
    
    if restore_button:
        try:
            # Validate UUID format
            UUID(version_id)
            
            result = _api_post(
                f"/api/v1/versions/restore/{version_id}",
                {
                    "confirm": True,
                    "reason": reason,
                }
            )
            
            st.success("✅ Version restored successfully!")
            st.info(f"**New version number**: {result['new_version_number']}")
            st.info(f"**Restored from version**: {result['restored_from_version']}")
            st.info(f"**Entity**: {result['entity_type']}/{result['entity_id']}")
            st.info(f"**Message**: {result['message']}")
            
            # Clear form
            st.session_state["restore_version_id"] = ""
            st.session_state["restore_reason"] = ""
            st.session_state["restore_confirm"] = False
            
        except ValueError:
            st.error("Invalid UUID format")
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                st.error("Version not found")
            elif e.response.status_code == 403:
                st.error("Admin access required")
            elif e.response.status_code == 400:
                st.error("Confirmation required")
            else:
                st.error(f"Failed to restore version: {e}")
        except Exception as e:
            st.error(f"Failed to restore version: {e}")


if __name__ == "__main__":
    render_provenance_ledger_page()
