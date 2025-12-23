"""Company settings (multi-tenancy) dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _headers() -> Dict[str, str]:
    # Companies endpoints use header-based auth (MVP): X-User-Id.
    uid = os.environ.get("DASHBOARD_USER_ID") or os.environ.get("X_USER_ID")
    try:
        # If Streamlit auth stores an id, prefer it.
        u = st.session_state.get("auth_user") or {}
        uid = u.get("id") or uid
    except Exception:
        pass
    return {"X-User-Id": str(uid)} if uid else {}


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_patch(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.patch(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_delete(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.delete(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json() if r.content else {"ok": True}


def render_company_settings_page() -> None:
    from scripts.dashboard.auth import require_auth, require_admin

    user = require_auth()
    require_admin(user)

    st.header("Company Settings")
    st.caption("Manage tenant company settings and users (superadmin only).")

    if not _headers():
        st.warning(
            "Missing `X-User-Id` context for Companies API. Set `DASHBOARD_USER_ID` env var to an admin user's UUID."
        )

    tab1, tab2 = st.tabs(["Company Info", "Users"])

    with tab1:
        st.subheader("Load company")
        company_id = st.text_input("Company UUID", value=st.session_state.get("company_settings_company_id", ""))
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Load"):
                st.session_state["company_settings_company_id"] = company_id
        with c2:
            if st.button("Create company"):
                st.session_state["company_settings_show_create"] = True

        if st.session_state.get("company_settings_show_create"):
            st.divider()
            st.subheader("Create company")
            name = st.text_input("Name", key="company_create_name")
            subdomain = st.text_input("Subdomain", key="company_create_subdomain")
            logo_url = st.text_input("Logo URL (optional)", key="company_create_logo")
            primary_color = st.text_input("Primary color (optional)", key="company_create_color", placeholder="#2E86AB")
            if st.button("Submit create", type="primary"):
                try:
                    out = _api_post(
                        "/api/companies",
                        {
                            "name": name,
                            "subdomain": subdomain,
                            "logo_url": logo_url or None,
                            "primary_color": primary_color or None,
                        },
                    )
                    st.success(f"Created company {out.get('id')}")
                    st.session_state["company_settings_company_id"] = out.get("id")
                    st.session_state["company_settings_show_create"] = False
                except Exception as e:  # noqa: BLE001
                    st.error(f"Create failed: {e}")

        st.divider()
        cid = st.session_state.get("company_settings_company_id")
        if cid:
            try:
                comp = _api_get(f"/api/companies/{cid}")
                st.json(comp)
                st.subheader("Edit company")
                name2 = st.text_input("Name", value=comp.get("name") or "", key="company_edit_name")
                logo2 = st.text_input("Logo URL", value=comp.get("logo_url") or "", key="company_edit_logo")
                color2 = st.text_input("Primary color", value=comp.get("primary_color") or "", key="company_edit_color")
                status2 = st.text_input("Status", value=comp.get("status") or "", key="company_edit_status")
                if st.button("Save changes"):
                    payload = {
                        "name": name2 or None,
                        "logo_url": logo2 or None,
                        "primary_color": color2 or None,
                        "status": status2 or None,
                    }
                    payload = {k: v for k, v in payload.items() if v is not None}
                    out = _api_patch(f"/api/companies/{cid}", payload)
                    st.success("Saved.")
                    st.json(out)
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load company: {e}")
        else:
            st.info("Enter a company UUID to load or create a new company.")

    with tab2:
        st.subheader("Users")
        cid = st.session_state.get("company_settings_company_id")
        if not cid:
            st.info("Load a company first.")
            return

        try:
            users = list(_api_get(f"/api/companies/{cid}/users") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load users: {e}")
            users = []

        if users:
            st.dataframe(pd.DataFrame(users), use_container_width=True, hide_index=True)
        else:
            st.info("No users found for this company.")

        st.divider()
        st.subheader("Invite user")
        email = st.text_input("Email", key="company_invite_email")
        role = st.selectbox("Company role", options=["member", "admin", "owner"], index=0, key="company_invite_role")
        if st.button("Invite", type="primary"):
            try:
                out = _api_post(f"/api/companies/{cid}/users/invite", {"email": email, "role": role}, timeout=60)
                st.success("Invited/assigned user.")
                st.json(out)
            except Exception as e:  # noqa: BLE001
                st.error(f"Invite failed: {e}")

        st.divider()
        st.subheader("Manage existing user")
        user_ids = [u.get("id") for u in users if u.get("id")]
        uid = st.selectbox("User", options=user_ids, index=0 if user_ids else None, key="company_manage_user")
        new_role = st.selectbox("New role", options=["member", "admin", "owner"], index=0, key="company_manage_role")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Change role"):
                try:
                    out = _api_patch(f"/api/companies/{cid}/users/{uid}", {"role": new_role}, timeout=60)
                    st.success("Role updated.")
                    st.json(out)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Role update failed: {e}")
        with c2:
            if st.button("Remove user", type="secondary"):
                try:
                    out = _api_delete(f"/api/companies/{cid}/users/{uid}", timeout=60)
                    st.success("User removed from company.")
                    st.json(out)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Remove failed: {e}")


__all__ = ["render_company_settings_page"]


