"""Lipidomics Spectral Matching dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List
from uuid import UUID

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_spectral_matching_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Spectral Matching")
    st.caption("Ingest spectral libraries and match lipid features using cosine similarity.")

    tab1, tab2 = st.tabs(["Library + Match", "Match Results"])

    with tab1:
        st.subheader("Library overview")
        libs: List[Dict[str, Any]] = []
        try:
            with httpx.Client(timeout=30) as client:
                r = client.get(f"{API_BASE}/api/spectral/libraries")
            if r.status_code == 200:
                libs = list(r.json() or [])
        except Exception:
            libs = []

        if libs:
            st.dataframe(pd.DataFrame(libs), use_container_width=True, hide_index=True)
        else:
            st.info("No spectral libraries found.")

        st.subheader("Ingest library (MGF)")
        mgf_path = st.text_input("MGF path", placeholder="/path/to/library.mgf")
        lib_name = st.text_input("Library name", placeholder="MyLibrary")
        lib_version = st.text_input("Version (optional)", placeholder="v1")
        if st.button("Ingest Library", type="primary"):
            payload = {"mgf_path": mgf_path, "name": lib_name, "version": lib_version or None}
            try:
                with httpx.Client(timeout=300) as client:
                    rr = client.post(f"{API_BASE}/api/spectral/libraries/ingest", json=payload)
                if rr.status_code >= 400:
                    st.error(f"Ingest failed ({rr.status_code}): {rr.text}")
                else:
                    st.success("Library ingested.")
            except Exception as e:  # noqa: BLE001
                st.error(f"Request failed: {e}")

        st.divider()
        st.subheader("Run match")
        feature_id = st.text_input("Feature UUID", placeholder="UUID of lipid feature to match")
        if st.button("Run Match"):
            try:
                fid = str(UUID(feature_id.strip()))
            except Exception:
                st.error("Invalid Feature UUID.")
                fid = ""
            if fid:
                try:
                    with httpx.Client(timeout=300) as client:
                        rr = client.post(f"{API_BASE}/api/spectral/match/{fid}")
                    if rr.status_code >= 400:
                        st.error(f"Match failed ({rr.status_code}): {rr.text}")
                    else:
                        st.success("Matching complete.")
                        st.session_state["spectral_last_feature_id"] = fid
                except Exception as e:  # noqa: BLE001
                    st.error(f"Request failed: {e}")

    with tab2:
        st.subheader("Match results")
        fid = st.text_input(
            "Feature UUID (results)",
            value=st.session_state.get("spectral_last_feature_id", ""),
            placeholder="UUID",
            key="spectral_results_feature",
        )
        if st.button("Refresh results"):
            st.session_state["spectral_last_feature_id"] = fid

        if not fid.strip():
            st.info("Enter a Feature UUID to view results.")
            return

        try:
            with httpx.Client(timeout=60) as client:
                rr = client.get(f"{API_BASE}/api/spectral/annotations/{fid.strip()}")
            if rr.status_code >= 400:
                st.error(f"Failed to load annotations ({rr.status_code}): {rr.text}")
                return
            rows = list(rr.json() or [])
        except Exception as e:  # noqa: BLE001
            st.error(f"Failed to load annotations: {e}")
            return

        if not rows:
            st.info("No annotations found (run matching first).")
            return

        df = pd.DataFrame(rows)
        cols = [
            "rank",
            "lipid_name",
            "lipid_class",
            "spectral_score",
            "mz_error_ppm",
            "is_confident",
            "is_ambiguous",
            "review_status",
        ]
        st.dataframe(df[[c for c in cols if c in df.columns]], use_container_width=True, hide_index=True)

        st.divider()
        st.subheader("Accept top hit (rename Feature)")
        top = rows[0]
        top_name = top.get("lipid_name") or ""
        if st.button(f'Accept "{top_name}"', disabled=not bool(top_name)):
            try:
                with httpx.Client(timeout=30) as client:
                    resp = client.patch(
                        f"{API_BASE}/api/v1/features/{fid.strip()}",
                        json={"name": top_name},
                    )
                if resp.status_code >= 400:
                    st.error(f"Rename failed ({resp.status_code}): {resp.text}")
                else:
                    st.success("Feature renamed.")
            except Exception as e:  # noqa: BLE001
                st.error(f"Rename request failed: {e}")

        st.caption("Review")
        for r in rows[:25]:
            ann_id = r.get("id")
            if not ann_id:
                continue
            with st.expander(f"Annotation {ann_id[:8]} rank={r.get('rank')} score={r.get('spectral_score')}"):
                st.json(r)
                c1, c2 = st.columns(2)
                with c1:
                    if st.button("Confirm", key=f"confirm_{ann_id}"):
                        with httpx.Client(timeout=30) as client:
                            resp = client.post(
                                f"{API_BASE}/api/spectral/annotations/{ann_id}/review",
                                json={"status": "confirmed"},
                            )
                        if resp.status_code == 200:
                            st.success("Confirmed.")
                        else:
                            st.error(f"Review failed ({resp.status_code}): {resp.text}")
                with c2:
                    if st.button("Reject", key=f"reject_{ann_id}"):
                        with httpx.Client(timeout=30) as client:
                            resp = client.post(
                                f"{API_BASE}/api/spectral/annotations/{ann_id}/review",
                                json={"status": "rejected"},
                            )
                        if resp.status_code == 200:
                            st.success("Rejected.")
                        else:
                            st.error(f"Review failed ({resp.status_code}): {resp.text}")


__all__ = ["render_spectral_matching_page"]


