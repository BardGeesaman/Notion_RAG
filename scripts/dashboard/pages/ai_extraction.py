"""AI Document Extraction dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post_multipart(path: str, files: List[tuple[str, tuple[str, bytes, str]]], *, timeout: int = 300) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", files=files)
    r.raise_for_status()
    return r.json()


def _api_get(path: str, *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _fetch_job(job_id: str) -> Dict[str, Any] | None:
    if not job_id:
        return None
    try:
        return _api_get(f"/api/extraction/jobs/{job_id}", timeout=60)
    except Exception:
        return None


def render_ai_extraction_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("AI Document Extraction")
    st.caption("Upload documents for parsing + structured entity extraction (LLM-backed).")

    tab1, tab2, tab3 = st.tabs(["Upload", "Jobs", "Results"])

    with tab1:
        st.subheader("Upload batch")
        uploads = st.file_uploader(
            "Select files",
            accept_multiple_files=True,
            type=["docx", "pptx", "xlsx", "csv", "pdf"],
        )

        if st.button("Start Extraction", type="primary", disabled=not uploads):
            multipart: List[tuple[str, tuple[str, bytes, str]]] = []
            for f in uploads or []:
                name = getattr(f, "name", "upload")
                content = f.getvalue()
                multipart.append(("files", (name, content, "application/octet-stream")))

            try:
                out = _api_post_multipart("/api/extraction/upload-batch", multipart, timeout=300)
                st.session_state["ai_extraction_job_id"] = out.get("job_id")
                st.session_state["ai_extraction_last_upload"] = out
            except Exception as e:  # noqa: BLE001
                st.error(f"Upload failed: {e}")

        last = st.session_state.get("ai_extraction_last_upload")
        job_id = st.session_state.get("ai_extraction_job_id") or (last.get("job_id") if isinstance(last, dict) else None)
        if last:
            st.json(last)

        if job_id:
            st.divider()
            st.subheader("Progress")
            col1, col2 = st.columns([3, 1])
            with col2:
                refresh = st.button("Refresh status")
            status = _fetch_job(job_id)
            if status:
                pct = float(status.get("progress_pct") or 0.0)
                st.progress(min(1.0, max(0.0, pct / 100.0)))
                st.caption(f"Job: {job_id} • {pct:.1f}%")
            else:
                st.info(f"Job: {job_id} (unable to fetch status)")

    with tab2:
        st.subheader("Recent jobs")
        skip = st.number_input("Skip", min_value=0, value=0, step=10)
        limit = st.number_input("Limit", min_value=1, max_value=200, value=25, step=5)
        if st.button("Load jobs"):
            try:
                listing = _api_get(f"/api/extraction/jobs?skip={int(skip)}&limit={int(limit)}", timeout=60)
                st.session_state["ai_extraction_jobs"] = listing
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load jobs: {e}")

        listing = st.session_state.get("ai_extraction_jobs") or {}
        jobs = listing.get("jobs") if isinstance(listing, dict) else None
        if isinstance(jobs, list) and jobs:
            st.caption(f"Total: {listing.get('total')} • Showing: {len(jobs)}")
            st.dataframe(jobs, use_container_width=True, hide_index=True)

            ids = [j.get("id") for j in jobs if isinstance(j, dict) and j.get("id")]
            chosen = st.selectbox("View job details", options=[""] + ids)
            if chosen:
                st.session_state["ai_extraction_job_id"] = chosen
                details = _fetch_job(chosen)
                if details:
                    st.json(details)
        else:
            st.info("Load jobs to view recent extraction runs.")

    with tab3:
        st.subheader("Results")
        # Prefer active job ID, else allow manual input
        active_job = st.session_state.get("ai_extraction_job_id") or ""
        job_id = st.text_input("Job ID", value=str(active_job))
        if st.button("Load results") and job_id:
            details = _fetch_job(job_id)
            if details:
                st.session_state["ai_extraction_job_details"] = details
            else:
                st.error("Failed to fetch job details.")

        details = st.session_state.get("ai_extraction_job_details")
        if isinstance(details, dict):
            st.caption(f"Progress: {details.get('progress_pct')}%")
            docs = details.get("documents") or []
            if isinstance(docs, list) and docs:
                for d in docs:
                    if not isinstance(d, dict):
                        continue
                    title = d.get("original_filename") or d.get("file_path") or d.get("id")
                    with st.expander(f"{title} • {d.get('status')}"):
                        st.json(d.get("extracted_entities") or {})
                        if d.get("error_log"):
                            st.caption(f"Error: {d.get('error_log')}")
            else:
                st.info("No documents found for this job.")


__all__ = ["render_ai_extraction_page"]


