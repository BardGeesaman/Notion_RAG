"""Docking Runs dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_docking_runs_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧪 Docking Runs")
    st.caption("Create and monitor Vina docking runs; review poses sorted by affinity.")

    # Load structures
    try:
        with httpx.Client(timeout=30) as client:
            sresp = client.get(f"{API_BASE}/api/structures", params={"limit": 200})
        sresp.raise_for_status()
        structures = list(sresp.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load structures: {e}")
        structures = []

    struct_map: Dict[str, str] = {}
    for s in structures:
        sid = s.get("id")
        if not sid:
            continue
        label = f"{sid[:8]} ({s.get('source')} {s.get('pdb_id') or s.get('alphafold_uniprot_id') or ''})"
        struct_map[label] = sid

    st.subheader("Create run")
    if not struct_map:
        st.info("No structures available. Fetch a structure first.")
        st.selectbox("Structure", ["(none)"])
        st.button("Start Docking Run", disabled=True)
    else:
        struct_label = st.selectbox("Structure", list(struct_map.keys()))
        structure_id = struct_map[struct_label]

        # Load pockets for selected structure
        pockets: List[Dict[str, Any]] = []
        try:
            with httpx.Client(timeout=30) as client:
                presp = client.get(f"{API_BASE}/api/pockets", params={"structure_id": structure_id})
            if presp.status_code == 200:
                pockets = list(presp.json() or [])
        except Exception:
            pockets = []

        pocket_map: Dict[str, str] = {"(none)": ""}
        for p in pockets:
            pid = p.get("id")
            if not pid:
                continue
            label = f"rank {p.get('pocket_rank')} score={p.get('score')} ({pid[:8]})"
            pocket_map[label] = pid

        pocket_label = st.selectbox("Binding site (optional)", list(pocket_map.keys()))
        binding_site_id = pocket_map[pocket_label] or None

        # Load compounds (string compound_id) via v1 API
        compounds: List[Dict[str, Any]] = []
        try:
            with httpx.Client(timeout=30) as client:
                cresp = client.get(f"{API_BASE}/api/v1/compounds/")
            if cresp.status_code == 200:
                compounds = list(cresp.json() or [])
        except Exception:
            compounds = []

        compound_choices = [c.get("compound_id") for c in compounds if c.get("compound_id")]
        selected_compounds = st.multiselect(
            "Compounds (multi-select)", compound_choices, default=compound_choices[:10]
        )
        st.caption("Note: compound selection uses compound_id strings; backend resolves to UUIDs.")

        if st.button("Start Docking Run", type="primary"):
            payload: Dict[str, Any] = {
                "structure_id": structure_id,
                "binding_site_id": binding_site_id,
                "compound_ids": selected_compounds,
            }
            try:
                with httpx.Client(timeout=60) as client:
                    rr = client.post(f"{API_BASE}/api/docking/runs", json=payload)
                if rr.status_code >= 400:
                    st.error(f"Create run failed ({rr.status_code}): {rr.text}")
                else:
                    st.success("Docking run started.")
            except Exception as e:  # noqa: BLE001
                st.error(f"Create run request failed: {e}")

    st.divider()
    st.subheader("Runs")

    status_filter = st.selectbox("Status filter", ["All", "pending", "running", "completed", "failed"])
    params: Dict[str, Any] = {"limit": 200}
    if status_filter != "All":
        params["status"] = status_filter

    try:
        with httpx.Client(timeout=30) as client:
            r = client.get(f"{API_BASE}/api/docking/runs", params=params)
        r.raise_for_status()
        runs = list(r.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load runs: {e}")
        return

    if not runs:
        st.info("No docking runs found.")
        return

    df = pd.DataFrame(
        [
            {
                "id": rr.get("id"),
                "status": rr.get("status"),
                "progress": float(rr.get("progress") or 0.0),
                "completed": rr.get("completed_compounds"),
                "total": rr.get("total_compounds"),
            }
            for rr in runs
        ]
    )
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.subheader("Run details")
    for rr in runs[:50]:
        rid = rr.get("id")
        if not rid:
            continue
        with st.expander(f"{rid} ({rr.get('status')})"):
            st.json(rr)
            st.progress(float(rr.get("progress") or 0.0))

            try:
                with httpx.Client(timeout=30) as client:
                    pr = client.get(f"{API_BASE}/api/docking/runs/{rid}/poses")
                if pr.status_code == 200:
                    poses = list(pr.json() or [])
                else:
                    poses = []
            except Exception:
                poses = []

            if poses:
                pdf = pd.DataFrame(
                    [
                        {
                            "pose_id": p.get("id"),
                            "compound_id": p.get("compound_id"),
                            "affinity": p.get("binding_affinity"),
                            "rmsd_lb": p.get("rmsd_lb"),
                            "rmsd_ub": p.get("rmsd_ub"),
                        }
                        for p in poses
                    ]
                )
                st.dataframe(pdf, use_container_width=True, hide_index=True)

                st.caption("Downloads")
                for p in poses[:25]:
                    pid = p.get("id")
                    if not pid:
                        continue
                    if st.button(f"Download pose {pid[:8]}", key=f"dl_pose_{pid}"):
                        with httpx.Client(timeout=30) as client:
                            resp = client.get(f"{API_BASE}/api/docking/poses/{pid}/download")
                        if resp.status_code == 200:
                            st.download_button(
                                label="Download PDBQT",
                                data=resp.content,
                                file_name=f"{pid}.pdbqt",
                                mime="chemical/x-pdbqt",
                                key=f"dl_pose_btn_{pid}",
                            )
                        else:
                            st.error(f"Download failed ({resp.status_code}): {resp.text}")
            else:
                st.info("No poses yet.")


__all__ = ["render_docking_runs_page"]


