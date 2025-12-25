"""Notebook Review Queue page (approval gates)."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st

from scripts.dashboard.components.notebook_card import render_notebook_card
from scripts.dashboard.components.review_card import render_review_badge


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _headers() -> Dict[str, str]:
    uid = os.environ.get("DASHBOARD_USER_ID") or os.environ.get("X_USER_ID")
    try:
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


def _api_put(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.put(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_review_queue_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Review Queue")
    st.caption("Approve/reject notebook dashboards with a signed review card.")

    if not _headers():
        st.warning("Missing `X-User-Id` context. Set `DASHBOARD_USER_ID` env var to a reviewer's UUID.")

    try:
        pending: List[Dict[str, Any]] = list(_api_get("/api/reviews/pending") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load pending reviews: {e}")
        return

    st.write(f"Pending reviews: {len(pending)}")
    if not pending:
        st.info("No pending reviews.")
        return

    for r in pending:
        rid = str(r.get("id") or "")
        nb_path = str(r.get("notebook_path") or "")
        with st.expander(f"{nb_path} • {rid}"):
            render_review_badge(
                {
                    "status": r.get("status"),
                    "reviewer_id": r.get("reviewer_id"),
                    "reviewed_at": r.get("reviewed_at"),
                    "signature": r.get("signature"),
                }
            )

            # Minimal preview card (actions only; diff is a follow-up)
            render_notebook_card(
                notebook_path=nb_path,
                title=nb_path,
                description=f"version_hash={r.get('version_hash')}",
                tags=[],
                voila_url=None,
                jupyter_url=None,
                source_url=None,
                preview_image=None,
            )

            comments = st.text_area("Comments", key=f"c_{rid}", height=100)
            c1, c2, c3 = st.columns(3)
            with c1:
                if st.button("Approve", type="primary", key=f"a_{rid}"):
                    try:
                        out = _api_put(f"/api/reviews/{rid}", {"status": "approved", "comments": comments})
                        st.success(f"Approved. Signature: {out.get('signature')}")
                        st.rerun()
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Approve failed: {e}")
            with c2:
                if st.button("Reject", key=f"r_{rid}"):
                    try:
                        _api_put(f"/api/reviews/{rid}", {"status": "rejected", "comments": comments})
                        st.success("Rejected.")
                        st.rerun()
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Reject failed: {e}")
            with c3:
                if st.button("Request changes", key=f"ch_{rid}"):
                    try:
                        _api_put(f"/api/reviews/{rid}", {"status": "changes_requested", "comments": comments})
                        st.success("Changes requested.")
                        st.rerun()
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Request changes failed: {e}")


__all__ = ["render_review_queue_page"]


