"""Authentication helpers for the dashboard."""

from __future__ import annotations

import streamlit as st
from amprenta_rag.auth.session import (
    get_current_user,
    clear_session,
    is_authenticated,
    update_last_activity,
    check_session_timeout,
    get_session_remaining,
    set_current_user,
)


def get_mock_user() -> dict:
    return {"id": "00000000-0000-0000-0000-000000000001", "username": "dev", "email": "dev@local", "role": "admin"}


def handle_session_timeout(auth_disabled: bool, timeout_minutes: int = 30) -> None:
    if auth_disabled or not is_authenticated():
        return
    if check_session_timeout(timeout_minutes=timeout_minutes):
        clear_session()
        st.warning("⏱️ Your session has expired. Please log in again.")
        from scripts.dashboard.pages.auth.login import render_login_page

        render_login_page()
        st.stop()
    else:
        remaining = get_session_remaining(timeout_minutes=timeout_minutes)
        if 0 < remaining < 5:
            st.warning(f"⏱️ Your session will expire in {remaining} minute{'s' if remaining != 1 else ''}. Please save your work.")


def check_authentication(auth_disabled: bool) -> dict | None:
    # Update last activity if authenticated
    if is_authenticated():
        update_last_activity()

    if not auth_disabled and not is_authenticated():
        from scripts.dashboard.pages.auth.login import render_login_page

        render_login_page()
        st.stop()

    if auth_disabled and not is_authenticated():
        set_current_user(get_mock_user())

    return get_current_user()

