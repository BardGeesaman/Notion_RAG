"""UI component: pin a notebook dashboard to a Program."""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


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


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_pin_dashboard_modal(*, program_id: UUID, key_prefix: str = "") -> Tuple[bool, Optional[str]]:
    """Render a popover 'modal' to pin a dashboard.

    Returns: (did_pin, error_message)
    """
    if not _headers():
        st.warning("Missing `X-User-Id` context. Set `DASHBOARD_USER_ID` env var to a user's UUID to pin dashboards.")
        return False, "missing_user_id"

    did_pin = False
    err: Optional[str] = None

    with st.popover("Pin Dashboard", use_container_width=False):
        st.subheader("Pin a notebook dashboard")

        try:
            notebooks: List[Dict[str, Any]] = list(_api_get("/api/notebooks") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load notebooks: {e}")
            return False, str(e)

        if not notebooks:
            st.info("No notebooks available to pin.")
            return False, None

        options = [(str(nb.get("title") or nb.get("path")), str(nb.get("path") or "")) for nb in notebooks]
        options = [(t, p) for (t, p) in options if p]
        if not options:
            st.info("No valid notebooks found.")
            return False, None

        labels = [f"{t} — {p}" for (t, p) in options]
        idx = 0
        selected = st.selectbox("Notebook", options=list(range(len(labels))), format_func=lambda i: labels[i], index=idx)
        notebook_path = options[int(selected)][1]

        display_name = st.text_input("Display name (optional)", key=f"{key_prefix}pin_display_name")
        config_text = st.text_area(
            "Voila params / config (JSON)",
            value="{}",
            height=120,
            key=f"{key_prefix}pin_config",
        )

        c1, c2 = st.columns(2)
        with c1:
            if st.button("Pin", type="primary", key=f"{key_prefix}pin_submit"):
                try:
                    cfg = json.loads(config_text) if (config_text or "").strip() else None
                    if cfg is not None and not isinstance(cfg, dict):
                        raise ValueError("config JSON must be an object")
                    _api_post(
                        "/api/dashboards/pin",
                        {
                            "notebook_path": notebook_path,
                            "program_id": str(program_id),
                            "display_name": display_name or None,
                            "config": cfg,
                        },
                    )
                    st.success("Pinned.")
                    did_pin = True
                except Exception as e:  # noqa: BLE001
                    err = str(e)
                    st.error(f"Pin failed: {e}")
        with c2:
            st.caption("Tip: set `DASHBOARD_USER_ID` for local dev auth.")

    return did_pin, err


__all__ = ["render_pin_dashboard_modal"]



