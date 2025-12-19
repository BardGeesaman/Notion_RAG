"""Streamlit session management for authentication."""
from typing import Optional
import time
import streamlit as st


def get_current_user() -> Optional[dict]:
    """Get the currently logged-in user from session state."""
    return st.session_state.get("user")


def set_current_user(user: dict) -> None:
    """Set the current user in session state."""
    st.session_state["user"] = user


def clear_session() -> None:
    """Clear the current user session (logout)."""
    st.session_state.pop("user", None)


def is_authenticated() -> bool:
    """Check if user is authenticated."""
    return get_current_user() is not None


def require_auth() -> bool:
    """Check authentication status. Returns False if not authenticated."""
    return is_authenticated()


def update_last_activity() -> None:
    """Update the last activity timestamp in session state."""
    st.session_state["last_activity"] = time.time()


def check_session_timeout(timeout_minutes: int = 30) -> bool:
    """
    Check if the session has timed out.

    Args:
        timeout_minutes: Session timeout in minutes

    Returns:
        True if session has timed out, False otherwise
    """
    if not is_authenticated():
        return False

    last_activity = st.session_state.get("last_activity")
    if last_activity is None:
        # Initialize if not set
        update_last_activity()
        return False

    elapsed_minutes = (time.time() - last_activity) / 60
    return elapsed_minutes >= timeout_minutes


def get_session_remaining(timeout_minutes: int = 30) -> int:
    """
    Get remaining session time in minutes.

    Args:
        timeout_minutes: Session timeout in minutes

    Returns:
        Minutes remaining until timeout (0 if timed out or not authenticated)
    """
    if not is_authenticated():
        return 0

    last_activity = st.session_state.get("last_activity")
    if last_activity is None:
        update_last_activity()
        return timeout_minutes

    elapsed_minutes = (time.time() - last_activity) / 60
    remaining = timeout_minutes - elapsed_minutes
    return max(0, int(remaining))
