"""Dashboard modules for the Streamlit application."""

from __future__ import annotations

from scripts.dashboard.db_session import db_session

# NOTE: Page render functions are imported lazily in run_dashboard.py
# to prevent cascading import failures. Don't add eager imports here.

__all__ = [
    "db_session",
]
