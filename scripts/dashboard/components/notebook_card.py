"""Reusable notebook card component for dashboard notebook listings."""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

import streamlit as st


JUPYTERHUB_URL = os.environ.get("JUPYTERHUB_URL", "http://localhost:8000").rstrip("/")
VOILA_BASE_URL = os.environ.get("VOILA_URL", "").rstrip("/") or JUPYTERHUB_URL


def render_notebook_card(
    *,
    notebook_path: str,
    title: str,
    description: str = "",
    tags: Optional[List[str]] = None,
    voila_url: Optional[str] = None,
    jupyter_url: Optional[str] = None,
    source_url: Optional[str] = None,
    preview_image: Optional[str] = None,
) -> None:
    """Render a single notebook card with action buttons."""

    nb = (notebook_path or "").strip()
    tags = tags or []

    st.subheader(title or Path(nb).stem or "Untitled")
    if description:
        st.caption(description)
    if tags:
        st.caption("Tags: " + ", ".join([str(t) for t in tags if t]))

    # Thumbnail (optional)
    if preview_image:
        try:
            p = Path(preview_image)
            if p.exists():
                st.image(str(p), use_container_width=True)
            else:
                st.caption("Preview: (placeholder)")
        except Exception:  # noqa: BLE001
            st.caption("Preview: (placeholder)")
    else:
        st.caption("Preview: (placeholder)")

    # Default links if not provided
    if jupyter_url is None and nb:
        jupyter_url = f"{JUPYTERHUB_URL}/notebooks/{nb}"
    if voila_url is None and nb:
        # Requirement: /voila/render/{notebook_path}
        voila_url = f"{VOILA_BASE_URL}/voila/render/{nb}"

    c1, c2, c3 = st.columns(3)
    with c1:
        if jupyter_url:
            st.link_button("Open in Jupyter", jupyter_url)
        else:
            st.button("Open in Jupyter", disabled=True)
    with c2:
        if voila_url:
            # Streamlit link buttons open in a new tab in most browsers.
            st.link_button("Open as Dashboard", voila_url)
        else:
            st.button("Open as Dashboard", disabled=True)
    with c3:
        if source_url:
            st.link_button("View Source", source_url)
        else:
            st.button("View Source", disabled=True)


__all__ = ["render_notebook_card"]


