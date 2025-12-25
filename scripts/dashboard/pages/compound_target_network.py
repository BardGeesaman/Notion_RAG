"""Compound-Target Network dashboard page."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Tuple

import httpx
import streamlit as st

from scripts.dashboard.components.cytoscape_network import DEFAULT_LAYOUTS, render_compound_target_network


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30, params: Dict[str, Any] | None = None) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}", params=params)
    r.raise_for_status()
    return r.json()


def _api_post(path: str, json_body: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=json_body)
    r.raise_for_status()
    return r.json()


def _log_slider_to_nm(v: float) -> float:
    return float(10 ** float(v))


def _format_nm(v: float) -> str:
    x = float(v)
    if x >= 1_000_000:
        return f"{x/1_000_000:.2f} mM"
    if x >= 1_000:
        return f"{x/1_000:.2f} µM"
    return f"{x:.1f} nM"


def _try_render_smiles_image(smiles: str) -> None:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except Exception:
        st.info("RDKit not installed; cannot render structure image.")
        return
    m = Chem.MolFromSmiles(str(smiles or ""))
    if m is None:
        st.error("Invalid SMILES.")
        return
    img = Draw.MolToImage(m, size=(420, 280))
    st.image(img, caption="Structure", use_container_width=False)


def render_compound_target_network_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Compound-Target Network")
    st.caption("Explore compound→target activity edges (GraphEdge relationship_type='activity_against').")

    tab_build, tab_vis, tab_details = st.tabs(["Build Network", "Visualize", "Details"])

    # Prefill from deep link / other pages
    prefill_compounds = st.session_state.pop("ctn_prefill_compound_ids", None)
    if prefill_compounds and "ctn_selected_compounds" not in st.session_state:
        st.session_state["ctn_selected_compounds"] = list(prefill_compounds)

    with tab_build:
        st.subheader("Build Network")

        # Load selector data (limit 200)
        if st.button("Refresh Lists", key="ctn_refresh_lists"):
            st.session_state.pop("ctn_compound_options", None)
            st.session_state.pop("ctn_target_options", None)

        comp_opts = st.session_state.get("ctn_compound_options")
        if not isinstance(comp_opts, list):
            try:
                comps = _api_get("/api/v1/compounds", timeout=30)
                comps = comps if isinstance(comps, list) else []
                comps = comps[:200]
                comp_opts = [{"id": c.get("id"), "compound_id": c.get("compound_id")} for c in comps if isinstance(c, dict)]
            except Exception as e:  # noqa: BLE001
                comp_opts = []
                st.error(f"Failed to load compounds: {e}")
            st.session_state["ctn_compound_options"] = comp_opts

        tgt_opts = st.session_state.get("ctn_target_options")
        if not isinstance(tgt_opts, list):
            try:
                feats = _api_get("/api/v1/features", params={"limit": 200, "feature_type": "gene"}, timeout=30)
                feats = feats if isinstance(feats, list) else []
                tgt_opts = [{"id": f.get("id"), "name": f.get("name")} for f in feats if isinstance(f, dict)]
            except Exception as e:  # noqa: BLE001
                tgt_opts = []
                st.error(f"Failed to load targets: {e}")
            st.session_state["ctn_target_options"] = tgt_opts

        comp_items: List[Tuple[str, str]] = []
        for c in comp_opts:
            cid = c.get("id")
            label = c.get("compound_id") or str(cid)[:8]
            if cid:
                comp_items.append((f"{label} ({str(cid)[:8]})", str(cid)))
        tgt_items: List[Tuple[str, str]] = []
        for t in tgt_opts:
            tid = t.get("id")
            nm = t.get("name") or str(tid)[:8]
            if tid:
                tgt_items.append((f"{nm} ({str(tid)[:8]})", str(tid)))

        comp_labels = [x[0] for x in comp_items]
        tgt_labels = [x[0] for x in tgt_items]
        comp_by_label = {lab: val for lab, val in comp_items}
        tgt_by_label = {lab: val for lab, val in tgt_items}

        default_comp_uuids = set(st.session_state.get("ctn_selected_compounds") or [])
        default_comp_labels = [lab for lab, uid in comp_items if uid in default_comp_uuids]

        sel_comp_labels = st.multiselect("Compounds (max 200)", options=comp_labels, default=default_comp_labels)
        sel_tgt_labels = st.multiselect("Targets (genes) (max 200)", options=tgt_labels, default=[])

        sel_comp_ids = [comp_by_label[x] for x in sel_comp_labels if x in comp_by_label]
        sel_tgt_ids = [tgt_by_label[x] for x in sel_tgt_labels if x in tgt_by_label]
        st.session_state["ctn_selected_compounds"] = sel_comp_ids
        st.session_state["ctn_selected_targets"] = sel_tgt_ids

        col_e1, col_e2 = st.columns(2)
        with col_e1:
            if st.button("Expand from Compounds", key="ctn_expand_compounds"):
                new_targets = set(sel_tgt_ids)
                for cid in sel_comp_ids:
                    try:
                        out = _api_get(f"/api/network/compound-target/expand/compound/{cid}", timeout=30)
                        tids = out.get("target_ids") if isinstance(out, dict) else []
                        for t in (tids or []):
                            new_targets.add(str(t))
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Expand failed for compound {cid[:8]}: {e}")
                st.session_state["ctn_selected_targets"] = list(new_targets)
                st.rerun()
        with col_e2:
            if st.button("Expand from Targets", key="ctn_expand_targets"):
                new_compounds = set(sel_comp_ids)
                for tid in sel_tgt_ids:
                    try:
                        out = _api_get(f"/api/network/compound-target/expand/target/{tid}", timeout=30)
                        cids = out.get("compound_ids") if isinstance(out, dict) else []
                        for c in (cids or []):
                            new_compounds.add(str(c))
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Expand failed for target {tid[:8]}: {e}")
                st.session_state["ctn_selected_compounds"] = list(new_compounds)
                st.rerun()

        st.markdown("### Filters")
        min_log, max_log = st.slider("IC50 range (log10 nM)", min_value=0.0, max_value=5.0, value=(0.0, 5.0), step=0.1)
        st.caption(f"Range: {_format_nm(_log_slider_to_nm(min_log))} → {_format_nm(_log_slider_to_nm(max_log))}")

        activity = st.selectbox("Activity type", options=["All", "IC50", "Ki", "Kd"], index=0)

        if st.button("Build Network", type="primary", key="ctn_build_network"):
            try:
                payload = {
                    "compound_ids": sel_comp_ids,
                    "target_ids": st.session_state.get("ctn_selected_targets") or [],
                    "filters": {
                        "ic50_range": {"min_nm": _log_slider_to_nm(min_log), "max_nm": _log_slider_to_nm(max_log)},
                        "activity_type": None if activity == "All" else activity,
                    },
                }
                out = _api_post("/api/network/compound-target", payload, timeout=60)
                st.session_state["ctn_network"] = out
            except Exception as e:  # noqa: BLE001
                st.session_state.pop("ctn_network", None)
                st.error(f"Network build failed: {e}")

        net = st.session_state.get("ctn_network")
        if isinstance(net, dict):
            st.caption(f"Nodes: {int((net.get('meta') or {}).get('node_count', 0) or 0)} · Edges: {int((net.get('meta') or {}).get('edge_count', 0) or 0)}")

    with tab_vis:
        st.subheader("Visualize")
        net = st.session_state.get("ctn_network")
        if not isinstance(net, dict) or not (net.get("nodes") and net.get("edges")):
            st.info("Build a network first (or no edges found).")
        else:
            layout = st.selectbox("Layout", options=DEFAULT_LAYOUTS, index=0, key="ctn_layout")
            st.caption("Legend: ellipse=compound · diamond=target · red=inhibitor · green=activator")
            st.caption(f"Node count: {len(net.get('nodes') or [])} (max 500)")
            selected = render_compound_target_network(
                nodes=net.get("nodes") or [],
                edges=net.get("edges") or [],
                layout=layout,
                options={},
                height=650,
                key="ctn_cy",
            )
            if isinstance(selected, dict):
                st.session_state["ctn_selected"] = selected

    with tab_details:
        st.subheader("Details")
        sel = st.session_state.get("ctn_selected")
        if not isinstance(sel, dict) or "kind" not in sel or "data" not in sel:
            st.info("Click a node or edge in the visualization to view details.")
            return

        kind = sel.get("kind")
        data = sel.get("data") if isinstance(sel.get("data"), dict) else {}

        if kind == "node":
            ntype = data.get("node_type") or data.get("entity_type")
            if ntype == "compound":
                st.markdown("### Compound")
                st.write(f"**Entity UUID:** `{data.get('entity_id')}`")
                st.write(f"**Compound ID:** `{data.get('compound_id')}`")
                # Best-effort fetch compound details for SMILES/MW if compound_id is present
                cid_str = data.get("compound_id")
                if cid_str:
                    try:
                        c = _api_get(f"/api/v1/compounds/{cid_str}", timeout=30)
                        if isinstance(c, dict):
                            smi = c.get("smiles")
                            st.write(f"**SMILES:** `{smi}`" if smi else "**SMILES:** -")
                            mw = c.get("molecular_weight")
                            st.write(f"**MW:** {mw}" if mw is not None else "**MW:** -")
                            if smi:
                                _try_render_smiles_image(str(smi))
                    except Exception as e:  # noqa: BLE001
                        st.caption(f"Compound fetch failed: {e}")
            else:
                st.markdown("### Target")
                st.write(f"**Entity UUID:** `{data.get('entity_id')}`")
                st.write(f"**Name:** `{data.get('target_name') or data.get('label')}`")
                st.write(f"**Feature type:** `{data.get('entity_type')}`")

        elif kind == "edge":
            st.markdown("### Edge")
            st.write(f"**Edge ID:** `{data.get('id')}`")
            st.write(f"**Best IC50 (nM):** `{data.get('best_ic50_nm')}`")
            st.write(f"**Activity type:** `{data.get('activity_type')}`")
            st.write(f"**Assay count:** `{data.get('assays_count')}`")
            assay_ids = data.get("assay_ids") or []
            if isinstance(assay_ids, list) and assay_ids:
                with st.expander(f"Assay IDs ({len(assay_ids)})", expanded=False):
                    st.code("\n".join([str(x) for x in assay_ids]))
            prov = data.get("provenance")
            if prov:
                with st.expander("Provenance", expanded=False):
                    st.json(prov)


__all__ = ["render_compound_target_network_page"]


