"""Executive Digests manager page (weekly Papermill + Voila/HTML digests)."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st

from amprenta_rag.database.models import Program
from scripts.dashboard.db_session import db_session


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


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_put(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.put(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def _api_delete(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout, headers=_headers()) as client:
        r = client.delete(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json() if r.content else {"ok": True}


def _cron_for_weekly(day_of_week: str, hour: int) -> str:
    return f"0 {int(hour)} * * {day_of_week}"


def render_digest_manager_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Executive Digests")
    st.caption("Schedule weekly notebook digests (Papermill execution + HTML output).")

    if not _headers():
        st.warning("Missing `X-User-Id` context. Set `DASHBOARD_USER_ID` env var to a user's UUID.")

    tab1, tab2 = st.tabs(["Schedules", "Create"])

    with tab1:
        st.subheader("Existing schedules")
        try:
            schedules: List[Dict[str, Any]] = list(_api_get("/api/digests") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load schedules: {e}")
            schedules = []

        if not schedules:
            st.info("No digest schedules yet.")
        else:
            for s in schedules:
                sid = str(s.get("id") or "")
                title = f"{s.get('notebook_path')} → program {s.get('program_id')}"
                with st.expander(title):
                    st.write(f"**Cron:** `{s.get('schedule_cron')}`")
                    st.write(f"**Recipients:** {', '.join(s.get('recipients') or [])}")
                    st.write(f"**Enabled:** {bool(s.get('enabled'))}")
                    st.write(f"**Last run:** {s.get('last_run_at')}")
                    st.write(f"**Last status:** {s.get('last_status')}")

                    c1, c2, c3 = st.columns(3)
                    with c1:
                        if st.button("Run now", key=f"run_{sid}"):
                            try:
                                out = _api_post(f"/api/digests/{sid}/run", {})
                                st.success(f"Run complete: {out.get('status')}")
                            except Exception as e:  # noqa: BLE001
                                st.error(f"Run failed: {e}")
                    with c2:
                        out_url = f"{API_BASE}/api/digests/{sid}/output"
                        st.link_button("View last output", out_url)
                    with c3:
                        if st.button("Delete", key=f"del_{sid}"):
                            try:
                                _api_delete(f"/api/digests/{sid}")
                                st.success("Deleted.")
                                st.rerun()
                            except Exception as e:  # noqa: BLE001
                                st.error(f"Delete failed: {e}")

    with tab2:
        st.subheader("Create schedule")
        with db_session() as db:
            progs = db.query(Program).order_by(Program.name.asc()).all()
        prog_map = {f"{p.name} ({p.id})": p.id for p in progs}
        if not prog_map:
            st.info("No programs available.")
            return
        prog_label = st.selectbox("Program", options=list(prog_map.keys()))
        program_id = prog_map[prog_label]

        try:
            notebooks: List[Dict[str, Any]] = list(_api_get("/api/notebooks") or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load notebooks: {e}")
            notebooks = []

        nb_options = {f"{n.get('title') or n.get('path')} — {n.get('path')}": str(n.get("path") or "") for n in notebooks}
        nb_options = {k: v for k, v in nb_options.items() if v}
        if not nb_options:
            st.info("No notebooks available in registry.")
            return

        nb_label = st.selectbox("Notebook", options=list(nb_options.keys()))
        notebook_path = nb_options[nb_label]

        st.caption("Weekly schedule helper")
        dcol, hcol = st.columns(2)
        with dcol:
            dow = st.selectbox("Day of week", options=["mon", "tue", "wed", "thu", "fri", "sat", "sun"], index=0)
        with hcol:
            hour = st.slider("Hour (0-23, UTC)", 0, 23, 9)
        cron = st.text_input("Cron expression", value=_cron_for_weekly(dow, hour))

        recipients_text = st.text_area("Recipients (one per line; email or identifier)", height=120, value="")
        recipients = [r.strip() for r in recipients_text.splitlines() if r.strip()]
        enabled = st.checkbox("Enabled", value=True)

        if st.button("Create", type="primary"):
            try:
                out = _api_post(
                    "/api/digests",
                    {
                        "program_id": str(program_id),
                        "notebook_path": notebook_path,
                        "schedule_cron": cron,
                        "recipients": recipients,
                        "enabled": enabled,
                    },
                )
                st.success(f"Created schedule {out.get('id')}")
                st.rerun()
            except Exception as e:  # noqa: BLE001
                st.error(f"Create failed: {e}")


__all__ = ["render_digest_manager_page"]


