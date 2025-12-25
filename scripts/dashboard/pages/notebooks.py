"""Notebooks page (template listing + Voila dashboard links)."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st

from scripts.dashboard.components.notebook_card import render_notebook_card
from scripts.dashboard.components.review_card import render_review_badge


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


@st.cache_data(ttl=30, show_spinner=False)
def _review_status(notebook_path: str) -> dict:
    try:
        return dict(_api_get(f"/api/reviews/notebook/{notebook_path}", timeout=10) or {})
    except Exception:
        return {}


def render_notebooks_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Notebooks")
    st.caption("Browse curated notebooks and open them read-only via Voila.")

    try:
        items: List[Dict[str, Any]] = list(_api_get("/api/notebooks") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load notebooks: {e}")
        return

    if not items:
        st.info("No notebooks available.")
        return

    cols = st.columns(3)
    for i, nb in enumerate(items):
        col = cols[i % 3]
        with col:
            render_review_badge(_review_status(str(nb.get("path") or "")))
            render_notebook_card(
                notebook_path=str(nb.get("path") or ""),
                title=str(nb.get("title") or nb.get("path") or "Untitled"),
                description=str(nb.get("description") or ""),
                tags=list(nb.get("tags") or []),
                voila_url=str(nb.get("voila_url") or ""),
                jupyter_url=str(nb.get("jupyter_url") or ""),
                source_url=str(nb.get("source_url") or ""),
                preview_image=str(nb.get("preview_image") or "") if nb.get("preview_image") else None,
            )


__all__ = ["render_notebooks_page"]


