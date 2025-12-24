"""External sync monitoring dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def render_sync_monitor_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("External Sync Monitor")
    st.caption("Run and monitor ChEMBL/PubChem sync jobs, and resolve conflicts.")

    tab_jobs, tab_conflicts, tab_stats = st.tabs(["Jobs", "Conflicts", "Stats"])

    with tab_jobs:
        st.subheader("Run sync")
        col1, col2 = st.columns([2, 1])
        with col1:
            source = st.selectbox("Source", options=["chembl", "pubchem"], index=0)
        with col2:
            sync_type = st.selectbox("Sync type", options=["incremental", "full"], index=0)

        if st.button("Run Sync", type="primary"):
            try:
                out = _api_post("/api/sync/run", {"source": source, "sync_type": sync_type}, timeout=30)
                st.session_state["sync_monitor_last_job"] = out
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to start sync: {e}")

        last = st.session_state.get("sync_monitor_last_job")
        if isinstance(last, dict) and last.get("job_id"):
            st.success(f"Started job: {last.get('job_id')}")

        st.divider()
        st.subheader("Recent jobs")
        colA, colB, colC = st.columns([1, 1, 2])
        with colA:
            skip = st.number_input("Skip", min_value=0, value=0, step=10)
        with colB:
            limit = st.number_input("Limit", min_value=1, max_value=200, value=25, step=5)
        with colC:
            filter_source = st.selectbox("Filter by source (optional)", options=["", "chembl", "pubchem"], index=0)

        if st.button("Load jobs"):
            try:
                qs = f"skip={int(skip)}&limit={int(limit)}"
                if filter_source:
                    qs += f"&source={filter_source}"
                listing = _api_get(f"/api/sync/jobs?{qs}", timeout=60)
                st.session_state["sync_monitor_jobs"] = listing
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load jobs: {e}")

        listing = st.session_state.get("sync_monitor_jobs") or {}
        jobs = listing.get("jobs") if isinstance(listing, dict) else None
        if isinstance(jobs, list) and jobs:
            st.caption(f"Total: {listing.get('total')} • Showing: {len(jobs)}")
            st.dataframe(jobs, use_container_width=True, hide_index=True)

            ids = [j.get("id") for j in jobs if isinstance(j, dict) and j.get("id")]
            chosen = st.selectbox("View job details", options=[""] + ids)
            if chosen:
                try:
                    details = _api_get(f"/api/sync/jobs/{chosen}", timeout=60)
                    st.json(details)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Failed to fetch job: {e}")
        else:
            st.info("Load jobs to view recent sync runs.")

    with tab_conflicts:
        st.subheader("Unresolved conflicts")
        col1, col2 = st.columns([1, 2])
        with col1:
            filter_source = st.selectbox("Source", options=["", "chembl", "pubchem"], index=0, key="conf_src")
        with col2:
            status = st.selectbox(
                "Resolution status", options=["pending", "auto_merged", "manual_override", "ignored"], index=0
            )

        if st.button("Load conflicts"):
            try:
                qs = f"skip=0&limit=100&resolution_status={status}"
                if filter_source:
                    qs += f"&source={filter_source}"
                out = _api_get(f"/api/sync/conflicts?{qs}", timeout=60)
                st.session_state["sync_monitor_conflicts"] = out
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load conflicts: {e}")

        conflicts_data = st.session_state.get("sync_monitor_conflicts") or {}
        conflicts = conflicts_data.get("conflicts") if isinstance(conflicts_data, dict) else None
        if isinstance(conflicts, list) and conflicts:
            st.caption(f"Total: {conflicts_data.get('total')} • Showing: {len(conflicts)}")
            st.dataframe(conflicts, use_container_width=True, hide_index=True)

            ids = [c.get("id") for c in conflicts if isinstance(c, dict) and c.get("id")]
            chosen = st.selectbox("Resolve conflict", options=[""] + ids)
            if chosen:
                colA, colB = st.columns([1, 1])
                with colA:
                    if st.button("Auto-merge"):
                        _api_post(f"/api/sync/conflicts/{chosen}/resolve", {"resolution": "auto_merged"}, timeout=30)
                        st.success("Resolved as auto_merged")
                with colB:
                    if st.button("Ignore"):
                        _api_post(f"/api/sync/conflicts/{chosen}/resolve", {"resolution": "ignored"}, timeout=30)
                        st.success("Resolved as ignored")
        else:
            st.info("Load conflicts to review unresolved sync issues.")

    with tab_stats:
        st.subheader("Stats (recent job history)")
        try:
            listing = _api_get("/api/sync/jobs?skip=0&limit=200", timeout=60)
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load jobs: {e}")
            return

        jobs = listing.get("jobs") if isinstance(listing, dict) else None
        if not isinstance(jobs, list) or not jobs:
            st.info("No sync jobs found.")
            return

        # Aggregate counts per source.
        agg: Dict[str, Dict[str, float]] = {}
        for j in jobs:
            if not isinstance(j, dict):
                continue
            src = (j.get("source") or "unknown").strip().lower()
            a = agg.setdefault(src, {"jobs": 0.0, "synced": 0.0, "updated": 0.0, "new": 0.0, "conflicts": 0.0})
            a["jobs"] += 1.0
            a["synced"] += float(j.get("records_synced") or 0)
            a["updated"] += float(j.get("records_updated") or 0)
            a["new"] += float(j.get("records_new") or 0)
            a["conflicts"] += float(j.get("conflicts_detected") or 0)

        rows: List[Dict[str, Any]] = []
        for src, a in sorted(agg.items()):
            rows.append({"source": src, **{k: int(v) for k, v in a.items()}})

        st.dataframe(rows, use_container_width=True, hide_index=True)


__all__ = ["render_sync_monitor_page"]


