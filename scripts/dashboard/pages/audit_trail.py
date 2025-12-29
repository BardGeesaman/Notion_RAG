"""Audit Trail Viewer page."""

from __future__ import annotations

import os
from typing import Any

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


def render_audit_trail_page() -> None:
    """Render the Audit Trail Viewer page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🔍 Audit Trail")
    st.caption("View audit history, verify data integrity, and manage electronic signatures.")
    
    tab1, tab2, tab3, tab4 = st.tabs(["Entity Audit", "Recent Activity", "Integrity Check", "Signatures"])
    
    with tab1:
        render_entity_audit_tab()
    
    with tab2:
        render_recent_activity_tab()
    
    with tab3:
        render_integrity_check_tab()
    
    with tab4:
        render_signatures_tab()


def render_entity_audit_tab() -> None:
    """Render entity audit history tab."""
    st.subheader("Entity Audit History")
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment", "Compound"],
        )
    
    with col2:
        entity_id = st.text_input("Entity ID", placeholder="UUID")
    
    if st.button("Fetch Audit Trail", type="primary", disabled=not entity_id):
        try:
            trail = _api_get(f"/api/v1/audit/{entity_type.lower()}/{entity_id}")
            st.session_state["audit_trail"] = trail
        except Exception as e:
            st.error(f"Failed to fetch audit trail: {e}")
    
    # Display trail
    trail = st.session_state.get("audit_trail")
    
    if trail:
        if trail:
            df = pd.DataFrame(trail)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.caption(f"Total entries: {len(trail)}")
        else:
            st.info("No audit history found for this entity")


def render_recent_activity_tab() -> None:
    """Render recent activity tab."""
    st.subheader("Recent Audit Activity")
    
    limit = st.slider("Entries to show", min_value=10, max_value=200, value=50, step=10)
    
    if st.button("Load Recent Activity"):
        try:
            recent = _api_get(f"/api/v1/audit/recent?limit={limit}")
            st.session_state["recent_activity"] = recent
        except Exception as e:
            st.error(f"Failed to load activity: {e}")
    
    # Display activity
    activity = st.session_state.get("recent_activity")
    
    if activity:
        if activity:
            df = pd.DataFrame(activity)
            
            # Select relevant columns
            display_cols = ["timestamp", "username", "action", "entity_type", "entity_id"]
            display_cols = [c for c in display_cols if c in df.columns]
            
            st.dataframe(df[display_cols], use_container_width=True, hide_index=True)
        else:
            st.info("No recent activity")


def render_integrity_check_tab() -> None:
    """Render integrity verification tab."""
    st.subheader("Data Integrity Verification")
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment", "Compound"],
            key="integrity_entity_type",
        )
    
    with col2:
        entity_id = st.text_input("Entity ID", placeholder="UUID", key="integrity_entity_id")
    
    # Current data input (simplified - would fetch from DB in production)
    st.markdown("**Current Data (JSON)**")
    current_data_str = st.text_area(
        "Paste current data JSON",
        height=100,
        placeholder='{"name": "test", "value": 123}',
        key="current_data",
    )
    
    if st.button("Verify Integrity", type="primary", disabled=not (entity_id and current_data_str)):
        try:
            import json
            current_data = json.loads(current_data_str)
            
            result = _api_post(
                f"/api/v1/audit/{entity_type.lower()}/{entity_id}/verify",
                {"current_data": current_data},
            )
            
            if result.get("verified"):
                st.success("✅ Data integrity verified!")
            else:
                st.error("❌ Data integrity check failed - checksum mismatch")
            
            st.info(result.get("message", ""))
            
        except json.JSONDecodeError:
            st.error("Invalid JSON format")
        except Exception as e:
            st.error(f"Verification failed: {e}")


def render_signatures_tab() -> None:
    """Render electronic signatures tab."""
    st.subheader("Electronic Signatures")
    
    # Sign action form
    st.markdown("### Create Signature")
    
    col1, col2 = st.columns(2)
    
    with col1:
        entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment", "Compound", "Protocol"],
            key="sig_entity_type",
        )
    
    with col2:
        entity_id = st.text_input("Entity ID", placeholder="UUID", key="sig_entity_id")
    
    action = st.selectbox(
        "Action",
        options=["Approve", "Reject", "Review", "Certify"],
        key="sig_action",
    )
    
    meaning = st.text_input(
        "Meaning",
        placeholder="e.g., I approve this compound for synthesis",
        key="sig_meaning",
    )
    
    password = st.text_input(
        "Password Confirmation",
        type="password",
        placeholder="Enter your password to sign",
        key="sig_password",
    )
    
    user_id = st.text_input(
        "User ID",
        placeholder="Your user UUID",
        key="sig_user_id",
        help="Get from your profile",
    )
    
    if st.button("✍️ Sign", type="primary", disabled=not all([entity_id, meaning, password, user_id])):
        try:
            result = _api_post(
                "/api/v1/signatures/sign",
                {
                    "user_id": user_id,
                    "action": action.lower(),
                    "entity_type": entity_type.lower(),
                    "entity_id": entity_id,
                    "meaning": meaning,
                    "password": password,
                },
            )
            
            st.success(f"✅ Signature created! ID: {result.get('signature_id')}")
            st.code(result.get("signature_hash", ""), language=None)
            
            # Refresh signature list
            st.session_state.pop("signature_list", None)
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 401:
                st.error("❌ Authentication failed - incorrect password")
            else:
                st.error(f"Failed to create signature: {e}")
        except Exception as e:
            st.error(f"Failed to create signature: {e}")
    
    # Signature history
    st.markdown("---")
    st.markdown("### Signature History")
    
    col1, col2, col3 = st.columns([2, 2, 1])
    
    with col1:
        hist_entity_type = st.selectbox(
            "Entity Type",
            options=["Dataset", "Experiment", "Compound"],
            key="hist_entity_type",
        )
    
    with col2:
        hist_entity_id = st.text_input("Entity ID", key="hist_entity_id")
    
    with col3:
        if st.button("Load Signatures"):
            if hist_entity_id:
                try:
                    sigs = _api_get(f"/api/v1/signatures/{hist_entity_type.lower()}/{hist_entity_id}")
                    st.session_state["signature_list"] = sigs
                except Exception as e:
                    st.error(f"Failed to load signatures: {e}")
    
    # Display signatures
    signature_list = st.session_state.get("signature_list")
    
    if signature_list:
        if signature_list:
            import pandas as pd
            
            df = pd.DataFrame(signature_list)
            
            # Truncate signature hash for display
            if "signature_hash" in df.columns:
                df["signature_hash"] = df["signature_hash"].str[:16] + "..."
            
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Verify signature
            if signature_list:
                sig_ids = [s.get("id") for s in signature_list if s.get("id")]
                if sig_ids:
                    selected_sig_id = st.selectbox("Select signature to verify", sig_ids)
                    
                    if st.button("🔒 Verify Signature"):
                        try:
                            verify_result = _api_get(f"/api/v1/signatures/{selected_sig_id}/verify")
                            
                            if verify_result.get("valid"):
                                st.success("✅ Signature is valid and untampered")
                            else:
                                st.error("❌ Signature verification failed - may be tampered")
                            
                        except Exception as e:
                            st.error(f"Verification failed: {e}")
        else:
            st.info("No signatures found for this entity")


if __name__ == "__main__":
    render_audit_trail_page()

