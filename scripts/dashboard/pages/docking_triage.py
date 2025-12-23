"""Docking Triage dashboard page (optional)."""

from __future__ import annotations

import os

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_docking_triage_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🎯 Docking Triage")
    st.caption("Filter docking hits by affinity and batch-download top poses.")

    try:
        with httpx.Client(timeout=30) as client:
            r = client.get(f"{API_BASE}/api/docking/runs", params={"limit": 200, "status": "completed"})
        r.raise_for_status()
        runs = list(r.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load completed runs: {e}")
        return

    if not runs:
        st.info("No completed runs found.")
        return

    run_map = {f"{rr.get('id')}": rr.get("id") for rr in runs if rr.get("id")}
    run_id = st.selectbox("Completed run", list(run_map.keys()))

    threshold = st.slider("Affinity threshold (kcal/mol)", min_value=-12.0, max_value=-4.0, value=-7.0, step=0.1)
    top_n = st.number_input("Top N to show/download", min_value=1, max_value=200, value=25, step=1)

    try:
        with httpx.Client(timeout=30) as client:
            pr = client.get(f"{API_BASE}/api/docking/runs/{run_id}/poses")
        pr.raise_for_status()
        poses = list(pr.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load poses: {e}")
        return

    hits = [p for p in poses if p.get("binding_affinity") is not None and float(p["binding_affinity"]) < float(threshold)]
    hits = hits[: int(top_n)]

    if not hits:
        st.info("No hits under threshold.")
        return

    df = pd.DataFrame(
        [
            {
                "pose_id": p.get("id"),
                "compound_id": p.get("compound_id"),
                "affinity": p.get("binding_affinity"),
                "rmsd_lb": p.get("rmsd_lb"),
                "rmsd_ub": p.get("rmsd_ub"),
            }
            for p in hits
        ]
    )
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.caption("Batch download")
    for p in hits:
        pid = p.get("id")
        if not pid:
            continue
        if st.button(f"Download {pid[:8]}", key=f"triage_dl_{pid}"):
            with httpx.Client(timeout=30) as client:
                resp = client.get(f"{API_BASE}/api/docking/poses/{pid}/download")
            if resp.status_code == 200:
                st.download_button(
                    label="Download PDBQT",
                    data=resp.content,
                    file_name=f"{pid}.pdbqt",
                    mime="chemical/x-pdbqt",
                    key=f"triage_dl_btn_{pid}",
                )
            else:
                st.error(f"Download failed ({resp.status_code}): {resp.text}")


__all__ = ["render_docking_triage_page"]


