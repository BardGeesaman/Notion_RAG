"""Protein Structures dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List
from uuid import UUID

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")
LIST_ENDPOINT = f"{API_BASE}/api/structures"
FETCH_ENDPOINT = f"{API_BASE}/api/structures/fetch"


def _get_structures(source: str | None = None, prep_status: str | None = None) -> List[Dict[str, Any]]:
    params: Dict[str, Any] = {"limit": 200}
    if source:
        params["source"] = source
    if prep_status:
        params["prep_status"] = prep_status
    with httpx.Client(timeout=30) as client:
        r = client.get(LIST_ENDPOINT, params=params)
        r.raise_for_status()
        return list(r.json() or [])


def render_protein_structures_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧬 Protein Structures")
    st.caption("Fetch, prepare, and manage protein structures (PDB / AlphaFold).")

    col1, col2 = st.columns(2)
    with col1:
        src_filter = st.selectbox("Filter source", ["All", "pdb", "alphafold"])
    with col2:
        prep_filter = st.selectbox("Filter prep status", ["All", "raw", "prepared", "failed"])

    st.subheader("Fetch structure")
    col_a, col_b, col_c = st.columns([1, 2, 2])
    with col_a:
        source = st.selectbox("Source", ["pdb", "alphafold"], key="struct_source")
    with col_b:
        pdb_id = st.text_input("PDB ID", placeholder="e.g., 1ABC", disabled=(source != "pdb"))
        uniprot_id = st.text_input("UniProt ID", placeholder="e.g., P12345", disabled=(source != "alphafold"))
    with col_c:
        feature_id = st.text_input("Feature UUID (optional)", placeholder="UUID")
        chain_id = st.text_input("Chain ID (optional)", placeholder="A")

    if st.button("Fetch", type="primary"):
        payload: Dict[str, Any] = {"source": source, "chain_id": chain_id or None}
        if source == "pdb":
            payload["pdb_id"] = pdb_id.strip()
        else:
            payload["uniprot_id"] = uniprot_id.strip()
        if feature_id.strip():
            try:
                payload["feature_id"] = str(UUID(feature_id.strip()))
            except Exception:
                st.error("Feature UUID invalid.")
                return
        try:
            with httpx.Client(timeout=60) as client:
                r = client.post(FETCH_ENDPOINT, json=payload)
            if r.status_code >= 400:
                st.error(f"Fetch failed ({r.status_code}): {r.text}")
                return
            st.success("Fetched structure.")
        except Exception as e:  # noqa: BLE001
            st.error(f"Request failed: {e}")

    st.divider()
    st.subheader("Existing structures")

    source_q = None if src_filter == "All" else src_filter
    prep_q = None if prep_filter == "All" else prep_filter
    try:
        rows = _get_structures(source=source_q, prep_status=prep_q)
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load structures: {e}")
        return

    if not rows:
        st.info("No structures found.")
        return

    df = pd.DataFrame(
        [
            {
                "id": r.get("id"),
                "source": r.get("source"),
                "pdb_id": r.get("pdb_id"),
                "alphafold_uniprot_id": r.get("alphafold_uniprot_id"),
                "prep_status": r.get("prep_status"),
                "files": len(r.get("files") or []),
            }
            for r in rows
        ]
    )
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.caption("Details")
    for r in rows[:50]:
        sid = r.get("id")
        with st.expander(f"{sid} ({r.get('source')}, {r.get('prep_status')})"):
            st.json(r)

            col_p1, col_p2, col_p3 = st.columns([1, 1, 2])
            with col_p1:
                do_prep = st.button("Prepare", key=f"prep_{sid}")
            with col_p2:
                prep_chain = st.text_input("Chain (optional)", key=f"prep_chain_{sid}")
            with col_p3:
                st.caption("Preparation uses pdbfixer/openmm on the backend.")

            if do_prep:
                try:
                    with httpx.Client(timeout=300) as client:
                        rr = client.post(f"{API_BASE}/api/structures/{sid}/prepare", json={"chain_id": prep_chain or None})
                    if rr.status_code >= 400:
                        st.error(f"Prepare failed ({rr.status_code}): {rr.text}")
                    else:
                        st.success("Prepared.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Prepare request failed: {e}")

            st.markdown("### Downloads")
            for ft in ["pdb", "prepared"]:
                try:
                    with httpx.Client(timeout=30) as client:
                        resp = client.get(f"{API_BASE}/api/structures/{sid}/download", params={"file_type": ft})
                    if resp.status_code == 200:
                        st.download_button(
                            f"Download {ft}",
                            data=resp.content,
                            file_name=f"{sid}_{ft}.pdb",
                            mime="chemical/x-pdb",
                            key=f"dl_{sid}_{ft}",
                        )
                except Exception:
                    continue


__all__ = ["render_protein_structures_page"]


