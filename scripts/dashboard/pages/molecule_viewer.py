"""Standalone 3D molecule viewer page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st

from scripts.dashboard.components.mol3d_viewer import render_conformers_3d, render_molecule_3d, render_overlay_3d


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 120) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def _parse_smiles_lines(text: str, max_n: int = 25) -> List[str]:
    lines = [ln.strip() for ln in (text or "").splitlines()]
    out = [ln for ln in lines if ln]
    return out[: int(max_n)]


def render_molecule_viewer_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Molecule Viewer")
    st.caption("Generate and visualize 3D conformers and overlays (best-effort, RDKit + py3Dmol optional).")

    tab1, tab2, tab3 = st.tabs(["Single Molecule", "Compare Molecules", "From Compound"])

    with tab1:
        st.subheader("Single Molecule")
        smi = st.text_input("SMILES", value="CCO", key="mv_smi")
        style = st.selectbox("Style", options=["stick", "sphere", "line"], index=0, key="mv_style")
        n_conf = st.slider("Conformer count", min_value=1, max_value=10, value=5, step=1, key="mv_nconf")

        if st.button("Generate", type="primary", key="mv_generate"):
            # Visual render (local)
            render_conformers_3d(smi, n_conformers=int(n_conf), style=style)
            # Energies + PDBs (API; used for download)
            try:
                out = _api_post(
                    "/api/viz3d/conformers",
                    {"smiles": smi, "n_conformers": int(n_conf), "optimize": True},
                    timeout=180,
                )
                st.session_state["mv_last_conformers"] = out
            except Exception as e:  # noqa: BLE001
                st.session_state.pop("mv_last_conformers", None)
                st.error(f"Failed to fetch conformers from API: {e}")

        cached = st.session_state.get("mv_last_conformers")
        if isinstance(cached, dict):
            pdbs = cached.get("pdb_strings") or []
            energies = cached.get("energies") or []
            if isinstance(energies, list) and energies:
                try:
                    df = pd.DataFrame({"conformer": list(range(len(energies))), "energy": energies})
                    st.dataframe(df, use_container_width=True, hide_index=True)
                except Exception:
                    st.write(energies)

            if isinstance(pdbs, list) and pdbs:
                idx = st.selectbox("Download conformer", options=list(range(len(pdbs))), index=0, key="mv_dl_idx")
                try:
                    pdb_str = str(pdbs[int(idx)])
                    st.download_button(
                        "Download PDB",
                        data=pdb_str.encode("utf-8"),
                        file_name="conformer.pdb",
                        mime="chemical/x-pdb",
                        key="mv_dl_pdb",
                    )
                except Exception:
                    pass

        st.divider()
        st.caption("Lowest-energy single view")
        render_molecule_3d(smi, style=style, width=420, height=420)

    with tab2:
        st.subheader("Compare Molecules")
        txt = st.text_area("SMILES list (one per line)", height=140, key="mv_compare_list")
        smiles_list = _parse_smiles_lines(txt, max_n=10)
        if smiles_list:
            ref_idx = st.selectbox("Reference", options=list(range(len(smiles_list))), index=0, key="mv_ref_idx")
            colors = []
            cols = st.columns(min(len(smiles_list), 5))
            for i in range(min(len(smiles_list), 5)):
                with cols[i]:
                    colors.append(st.color_picker(f"Mol {i}", value=["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"][i], key=f"mv_c_{i}"))
            if st.button("Overlay", type="primary", key="mv_overlay"):
                render_overlay_3d(smiles_list, colors=colors, reference_idx=int(ref_idx))
        else:
            st.info("Enter 2+ SMILES to overlay.")

    with tab3:
        st.subheader("From Compound")
        cid = st.text_input("Compound ID", placeholder="e.g., CMPD-0001", key="mv_compound_id")
        if st.button("Load Compound", key="mv_load_compound", type="secondary"):
            try:
                out = _api_get(f"/api/v1/compounds/{cid}", timeout=30)
                if isinstance(out, dict) and out.get("smiles"):
                    st.session_state["mv_compound_smiles"] = str(out["smiles"])
                else:
                    st.session_state.pop("mv_compound_smiles", None)
                    st.error("Compound found but has no SMILES.")
            except Exception as e:  # noqa: BLE001
                st.session_state.pop("mv_compound_smiles", None)
                st.error(f"Failed to load compound: {e}")

        csmi = st.session_state.get("mv_compound_smiles")
        if isinstance(csmi, str) and csmi.strip():
            st.write(f"**SMILES:** `{csmi}`")
            style2 = st.selectbox("Style", options=["stick", "sphere", "line"], index=0, key="mv_c_style")
            n_conf2 = st.slider("Conformer count", min_value=1, max_value=10, value=5, step=1, key="mv_c_nconf")
            if st.button("Generate 3D", type="primary", key="mv_c_generate"):
                render_conformers_3d(csmi, n_conformers=int(n_conf2), style=style2)
            render_molecule_3d(csmi, style=style2, width=420, height=420)


__all__ = ["render_molecule_viewer_page"]


