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
    st.caption("Filter docking hits by affinity, run pose QC (PLIP), and batch-download poses.")

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

    # Fetch QC summaries (best-effort) to enrich the table
    if "pose_qc_cache" not in st.session_state:
        st.session_state["pose_qc_cache"] = {}
    qc_cache = st.session_state["pose_qc_cache"]

    with httpx.Client(timeout=30) as client:
        for p in hits:
            pid = p.get("id")
            if not pid or pid in qc_cache:
                continue
            resp = client.get(f"{API_BASE}/api/poses/{pid}/quality")
            if resp.status_code == 200:
                qc_cache[pid] = resp.json()
            else:
                qc_cache[pid] = None

    df["LE"] = [((qc_cache.get(p.get("id")) or {}).get("ligand_efficiency")) for p in hits]
    df["Interactions"] = [((qc_cache.get(p.get("id")) or {}).get("total_interactions")) for p in hits]
    df["QC Status"] = ["analyzed" if qc_cache.get(p.get("id")) else "not analyzed" for p in hits]

    st.dataframe(df, use_container_width=True, hide_index=True)

    st.subheader("Pose QC")
    st.caption("Analyze uses PLIP on the backend; results are cached per pose.")

    if "pose_interactions_cache" not in st.session_state:
        st.session_state["pose_interactions_cache"] = {}
    ix_cache = st.session_state["pose_interactions_cache"]

    for p in hits:
        pid = p.get("id")
        if not pid:
            continue
        aff = p.get("binding_affinity")
        with st.expander(f"Pose {pid[:8]}  affinity={aff}"):
            cols = st.columns([1, 1, 2])
            with cols[0]:
                do_analyze = st.button("Analyze", key=f"analyze_{pid}")
            with cols[1]:
                st.caption(f"QC: {('analyzed' if qc_cache.get(pid) else 'not analyzed')}")
            with cols[2]:
                st.caption("Runs PLIP + clash check + ligand efficiency.")

            if do_analyze:
                try:
                    with httpx.Client(timeout=300) as client:
                        rr = client.post(f"{API_BASE}/api/poses/{pid}/analyze")
                    if rr.status_code >= 400:
                        st.error(f"Analyze failed ({rr.status_code}): {rr.text}")
                    else:
                        qc_cache[pid] = rr.json()
                        ix_cache.pop(pid, None)
                        st.success("Analyzed.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Analyze request failed: {e}")

            qc = qc_cache.get(pid)
            if qc:
                st.json(qc)
                if st.button("Load interactions", key=f"load_ix_{pid}"):
                    with httpx.Client(timeout=60) as client:
                        ir = client.get(f"{API_BASE}/api/poses/{pid}/interactions")
                    if ir.status_code == 200:
                        ix_cache[pid] = list(ir.json() or [])
                    else:
                        st.error(f"Failed to load interactions ({ir.status_code}): {ir.text}")

                ix = ix_cache.get(pid)
                if ix:
                    st.caption("Interactions")
                    st.dataframe(pd.DataFrame(ix), use_container_width=True, hide_index=True)
            else:
                st.info("No QC data yet. Click Analyze to run Pose QC.")

    st.divider()
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


