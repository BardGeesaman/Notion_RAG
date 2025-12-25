"""Review status badge for notebook cards."""

from __future__ import annotations

from typing import Any, Dict

import streamlit as st


def render_review_badge(status_payload: Dict[str, Any]) -> None:
    status = (status_payload or {}).get("status")
    reviewer = (status_payload or {}).get("reviewer_id")
    reviewed_at = (status_payload or {}).get("reviewed_at")
    signature = (status_payload or {}).get("signature")

    if not status:
        st.caption("Review: (none)")
        return

    if status == "approved":
        st.success(f"Approved — sig {str(signature)[:12]}…")
    elif status == "pending":
        st.warning("Pending review")
    elif status in ("rejected", "changes_requested"):
        st.error(f"{status.replace('_', ' ').title()}")
    else:
        st.caption(f"Review: {status}")

    hover = []
    if reviewer:
        hover.append(f"reviewer={reviewer}")
    if reviewed_at:
        hover.append(f"at={reviewed_at}")
    if signature and status == "approved":
        hover.append(f"signature={signature}")
    if hover:
        st.caption(" • ".join(hover))


__all__ = ["render_review_badge"]


