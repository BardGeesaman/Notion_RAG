"""Dynamic routing for dashboard pages."""

from __future__ import annotations

import importlib
import streamlit as st

from .config import PAGE_REGISTRY


def route_to_page(page_name: str) -> None:
    entry = PAGE_REGISTRY.get(page_name)
    if not entry:
        st.error(f"Unknown page: {page_name}")
        return
    module_path, func_name = entry
    try:
        module = importlib.import_module(module_path)
        render_func = getattr(module, func_name, None)
        if not render_func:
            raise AttributeError(f"{func_name} not found in {module_path}")
        render_func()
    except ImportError as e:
        st.error(f"❌ Error loading page: {page_name}")
        st.error(f"Import error: {str(e)}")
        st.info("This page may depend on modules that are currently unavailable. Try another page.")
        with st.expander("Show full error details"):
            st.exception(e)
    except Exception as e:
        st.error(f"❌ Error rendering page: {page_name}")
        st.exception(e)

