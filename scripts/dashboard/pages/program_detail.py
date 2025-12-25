"""Program detail helpers for Streamlit dashboard."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional
from uuid import UUID

import httpx
import streamlit as st

from scripts.dashboard.components.notebook_card import render_notebook_card
from scripts.dashboard.components.pin_dashboard_modal import render_pin_dashboard_modal


API_BASE = os.environ.get("API_URL", "http://localhost:8000")
JUPYTERHUB_URL = os.environ.get("JUPYTERHUB_URL", "http://localhost:8000").rstrip("/")
VOILA_BASE_URL = os.environ.get("VOILA_URL", "").rstrip("/") or JUPYTERHUB_URL


def _headers() -> Dict[str, str]:
    uid = os.environ.get("DASHBOARD_USER_ID") or os.environ.get("X_USER_ID")
    try:
        u = st.session_state.get("auth_user") or {}
        uid = u.get("id") or uid
    except Exception:
        pass
    return {"X-User-Id": str(uid)} if uid else {}


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_delete(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.delete(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json() if r.content else {"ok": True}


def render_pinned_dashboards_section(*, program_id: UUID, key_prefix: str = "") -> None:
    """Render the Pinned Dashboards section inside a program view."""
    st.subheader("Pinned Dashboards")

    did_pin, _ = render_pin_dashboard_modal(program_id=program_id, key_prefix=f"{key_prefix}pin_")
    if did_pin:
        st.rerun()

    if not _headers():
        st.caption("Set `DASHBOARD_USER_ID` to enable pin/unpin actions.")
        return

    try:
        pins: List[Dict[str, Any]] = list(_api_get(f"/api/dashboards/program/{program_id}") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load pinned dashboards: {e}")
        return

    if not pins:
        st.info("No pinned dashboards for this program yet.")
        return

    cols = st.columns(3)
    for i, pin in enumerate(pins):
        col = cols[i % 3]
        with col:
            nb_path = str(pin.get("notebook_path") or "")
            title = str(pin.get("display_name") or nb_path or "Notebook")

            # Best-effort metadata enrichment
            desc: str = ""
            tags: Optional[List[str]] = None
            source_url: Optional[str] = None
            try:
                meta = _api_get(f"/api/notebooks/{nb_path}/metadata")
                if isinstance(meta, dict):
                    desc = str(meta.get("description") or "")
                    tags = list(meta.get("tags") or [])
                    src = meta.get("source_url")
                    if isinstance(src, str) and src.strip():
                        source_url = f"{API_BASE}{src}" if src.startswith("/") else src
                    if not pin.get("display_name"):
                        title = str(meta.get("title") or title)
            except Exception:
                pass

            render_notebook_card(
                notebook_path=nb_path,
                title=title,
                description=desc,
                tags=tags or [],
                jupyter_url=f"{JUPYTERHUB_URL}/notebooks/{nb_path}" if nb_path else None,
                voila_url=f"{VOILA_BASE_URL}/voila/render/{nb_path}" if nb_path else None,
                source_url=source_url,
                preview_image=None,
            )

            # Replace placeholder "View Source" with Unpin control for pinned cards.
            pin_id = str(pin.get("id") or "")
            if pin_id:
                if st.button("Unpin", key=f"{key_prefix}unpin_{pin_id}"):
                    try:
                        _api_delete(f"/api/dashboards/pin/{pin_id}")
                        st.success("Unpinned.")
                        st.rerun()
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Unpin failed: {e}")


__all__ = ["render_pinned_dashboards_section"]


