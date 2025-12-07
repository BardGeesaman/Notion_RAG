"""
Page modules for the Streamlit dashboard.

NOTE: This package uses LAZY IMPORTS to prevent cascading failures.
Each page is imported only when selected by the user in run_dashboard.py.

DO NOT add eager imports here - they cause the entire dashboard to crash
if any single page has an import error.
"""

from __future__ import annotations

# No eager imports - all pages loaded lazily in run_dashboard.py

__all__ = []
