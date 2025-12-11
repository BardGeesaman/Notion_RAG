"""Streamlit session management for authentication."""
from typing import Optional
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
