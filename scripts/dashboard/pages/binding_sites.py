"""Binding site (pocket) detection dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def render_binding_sites_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🕳️ Binding Sites")
    st.caption("Detect and review binding pockets (fpocket).")

    # Load structures via API
    structures: List[Dict[str, Any]] = []
    try:
        with httpx.Client(timeout=30) as client:
            r = client.get(f"{API_BASE}/api/structures", params={"limit": 200})
        r.raise_for_status()
        structures = list(r.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load structures: {e}")

    structure_id = None
    if structures:
        labels = {}
        for s in structures:
            sid = s.get("id")
            if not sid:
                continue
            label = f"{sid[:8]} ({s.get('source')} {s.get('pdb_id') or s.get('alphafold_uniprot_id') or ''})"
            labels[label] = sid
        if labels:
            selected_label = st.selectbox("Structure", list(labels.keys()))
            structure_id = labels[selected_label]
    else:
        st.selectbox("Structure", ["(no structures available)"])

    if st.button("Detect Pockets", type="primary"):
        if not structure_id:
            st.error("No structure selected.")
            return
        try:
            with httpx.Client(timeout=300) as client:
                rr = client.post(f"{API_BASE}/api/pockets/structures/{structure_id}/detect-pockets")
            if rr.status_code >= 400:
                st.error(f"Detection failed ({rr.status_code}): {rr.text}")
            else:
                st.success("Pocket detection completed.")
        except Exception as e:  # noqa: BLE001
            st.error(f"Detection request failed: {e}")

    st.divider()
    st.subheader("Pockets")

    if not structure_id:
        st.info("Select a structure to view pockets.")
        return

    try:
        with httpx.Client(timeout=30) as client:
            pr = client.get(f"{API_BASE}/api/pockets", params={"structure_id": structure_id})
        pr.raise_for_status()
        pockets = list(pr.json() or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load pockets: {e}")
        return

    if not pockets:
        st.info("No pockets detected yet.")
        return

    df = pd.DataFrame(
        [
            {
                "id": p.get("id"),
                "rank": p.get("pocket_rank"),
                "score": p.get("score"),
                "volume": p.get("volume"),
                "center": f"({p.get('center_x'):.2f}, {p.get('center_y'):.2f}, {p.get('center_z'):.2f})"
                if p.get("center_x") is not None
                else None,
                "method": p.get("detection_method"),
            }
            for p in pockets
        ]
    )
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.caption("Downloads")
    for p in pockets:
        pid = p.get("id")
        if not pid:
            continue
        if st.button(f"Download pocket {p.get('pocket_rank')}", key=f"dl_{pid}"):
            try:
                with httpx.Client(timeout=30) as client:
                    resp = client.get(f"{API_BASE}/api/pockets/{pid}/download")
                if resp.status_code >= 400:
                    st.error(f"Download failed ({resp.status_code}): {resp.text}")
                else:
                    st.download_button(
                        label="Download PDB",
                        data=resp.content,
                        file_name=f"pocket_{p.get('pocket_rank')}.pdb",
                        mime="chemical/x-pdb",
                        key=f"dl_btn_{pid}",
                    )
            except Exception as e:  # noqa: BLE001
                st.error(f"Download error: {e}")


__all__ = ["render_binding_sites_page"]


