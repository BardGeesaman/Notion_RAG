"""Pathway Map Viewer dashboard page (KEGG-only MVP)."""

from __future__ import annotations

import base64
import os
from typing import Any, Dict, List, Optional, Tuple

import httpx
import pandas as pd
import streamlit as st

from scripts.dashboard.components.cytoscape_pathway import render_color_legend, render_pathway_map


def _api_base() -> str:
    return str(st.session_state.get("API_URL") or os.environ.get("API_URL") or "http://localhost:8000").rstrip("/")


def _api_get(path: str, *, params: Optional[Dict[str, Any]] = None, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{_api_base()}{path}", params=params)
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 60) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{_api_base()}{path}", json=payload)
    r.raise_for_status()
    return r.json()


@st.cache_data(ttl=60)
def _list_datasets(limit: int = 200) -> List[Dict[str, Any]]:
    try:
        out = _api_get("/api/v1/datasets", params={"limit": int(limit)}, timeout=30)
        return out if isinstance(out, list) else []
    except Exception:
        return []


def _dataset_options() -> List[Tuple[str, str]]:
    dsets = _list_datasets(200)
    opts = []
    for d in dsets:
        if not isinstance(d, dict):
            continue
        did = d.get("id")
        name = d.get("name") or str(did)[:8]
        om = d.get("omics_type") or ""
        if did:
            opts.append((str(did), f"{name} {f'({om})' if om else ''}".strip()))
    return opts


def _ensure_structure_loaded(pathway_id: str) -> Optional[Dict[str, Any]]:
    pid = (pathway_id or "").strip()
    if not pid:
        return None
    try:
        stt = _api_get(f"/api/pathway-maps/structure/{pid}", timeout=60)
        if isinstance(stt, dict):
            st.session_state["pmv_structure"] = stt
            st.session_state["pmv_pathway_id"] = pid
            st.session_state["pmv_pathway_name"] = stt.get("name") or pid
            return stt
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load pathway structure: {e}")
    return None


def render_pathway_map_viewer_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Pathway Map Viewer")
    st.caption("KEGG KGML-based pathway visualization with optional gene expression overlays.")

    # Defaults
    st.session_state.setdefault("pmv_pathway_id", "")
    st.session_state.setdefault("pmv_overlays", [])
    st.session_state.setdefault("pmv_expression", {})
    st.session_state.setdefault("pmv_vmin", -2.0)
    st.session_state.setdefault("pmv_vmax", 2.0)
    st.session_state.setdefault("pmv_colormap", "RdBu_r")

    tab1, tab2, tab3, tab4 = st.tabs(["Browse Pathways", "Pathway Map", "Data Overlay", "Multi-Pathway"])

    with tab1:
        st.subheader("Browse Pathways")

        q = st.text_input("Search KEGG pathways", value=st.session_state.get("pmv_search_query", ""), key="pmv_search_query")
        if st.button("Search", type="primary", key="pmv_search_btn"):
            try:
                res = _api_get("/api/pathway-maps/search", params={"query": q}, timeout=30)
                st.session_state["pmv_search_results"] = res if isinstance(res, list) else []
            except Exception as e:  # noqa: BLE001
                st.session_state["pmv_search_results"] = []
                st.error(f"Search failed: {e}")

        search_results = st.session_state.get("pmv_search_results") or []
        if isinstance(search_results, list) and search_results:
            df = pd.DataFrame(search_results)
            for col in ("pathway_id", "name", "organism", "gene_count"):
                if col not in df.columns:
                    df[col] = None
            st.dataframe(df[["pathway_id", "name", "organism", "gene_count"]], hide_index=True, use_container_width=True)
            pid = st.selectbox("Select pathway from search results", options=[r.get("pathway_id") for r in search_results if isinstance(r, dict) and r.get("pathway_id")], key="pmv_pid_from_search")
            if st.button("Load Selected Pathway", key="pmv_load_from_search"):
                _ensure_structure_loaded(str(pid))

        st.divider()
        st.markdown("### Enriched Pathways (from dataset)")
        ds_opts = _dataset_options()
        if ds_opts:
            ds_id = st.selectbox("Dataset", options=[x[0] for x in ds_opts], format_func=lambda x: dict(ds_opts).get(x, x), key="pmv_enrich_ds")
            if st.button("Get Enriched", key="pmv_get_enriched"):
                try:
                    res = _api_get(f"/api/pathway-maps/enriched/{ds_id}", timeout=120)
                    st.session_state["pmv_enriched_results"] = res if isinstance(res, list) else []
                except Exception as e:  # noqa: BLE001
                    st.session_state["pmv_enriched_results"] = []
                    st.error(f"Enriched fetch failed: {e}")
        else:
            st.info("No datasets available.")

        enr = st.session_state.get("pmv_enriched_results") or []
        if isinstance(enr, list) and enr:
            df = pd.DataFrame(enr)
            for col in ("pathway_id", "name", "organism", "gene_count"):
                if col not in df.columns:
                    df[col] = None
            st.dataframe(df[["pathway_id", "name", "organism", "gene_count"]], hide_index=True, use_container_width=True)
            pid = st.selectbox("Select pathway from enriched results", options=[r.get("pathway_id") for r in enr if isinstance(r, dict) and r.get("pathway_id")], key="pmv_pid_from_enriched")
            if st.button("Load Enriched Pathway", key="pmv_load_from_enriched"):
                _ensure_structure_loaded(str(pid))

    with tab2:
        st.subheader("Pathway Map")
        pid = st.session_state.get("pmv_pathway_id") or ""
        pname = st.session_state.get("pmv_pathway_name") or pid or "(none)"
        st.markdown(f"**Selected pathway:** `{pname}`")

        layout_choice = st.selectbox("Layout", options=["KEGG Layout", "Force-Directed"], index=0, key="pmv_layout")
        layout = "preset" if layout_choice == "KEGG Layout" else "cose"

        stt = st.session_state.get("pmv_structure")
        if not isinstance(stt, dict) or not stt.get("nodes"):
            st.info("No pathway loaded yet. Use the Browse tab to search and select a pathway.")
        else:
            nodes = stt.get("nodes") or []
            edges = stt.get("edges") or []
            overlays = st.session_state.get("pmv_overlays") or []

            st.caption("Legend: gene=rectangle · compound=ellipse · sub-pathway=hexagon")
            val = render_pathway_map(nodes=nodes, edges=edges, overlays=overlays if isinstance(overlays, list) else None, layout=layout, height=700, key="pmv_cy")
            if isinstance(val, dict) and val.get("kind") == "png" and val.get("data_url"):
                st.session_state["pmv_last_png"] = val["data_url"]

            png = st.session_state.get("pmv_last_png")
            if isinstance(png, str) and png.startswith("data:image/png;base64,"):
                b64 = png.split(",", 1)[1]
                try:
                    raw = base64.b64decode(b64.encode("utf-8"))
                    st.download_button("Download PNG", data=raw, file_name=f"{pid or 'pathway'}.png", mime="image/png")
                except Exception:
                    st.caption("PNG export captured but could not decode.")
            else:
                st.caption("PNG export available via the 'Export PNG' button above the map.")

    with tab3:
        st.subheader("Data Overlay")
        pid = st.session_state.get("pmv_pathway_id") or ""
        if not pid:
            st.info("Select a pathway first (Browse tab).")
        else:
            ds_opts = _dataset_options()
            if not ds_opts:
                st.info("No datasets available for overlays.")
            else:
                ds_id = st.selectbox("Dataset", options=[x[0] for x in ds_opts], format_func=lambda x: dict(ds_opts).get(x, x), key="pmv_expr_ds")

                if st.button("Load Expression", key="pmv_load_expr"):
                    try:
                        out = _api_get(f"/api/pathway-maps/expression/{ds_id}", timeout=60)
                        expr = out.get("gene_expression") if isinstance(out, dict) else {}
                        st.session_state["pmv_expression"] = expr if isinstance(expr, dict) else {}
                    except Exception as e:  # noqa: BLE001
                        st.session_state["pmv_expression"] = {}
                        st.error(f"Expression load failed: {e}")

                expr = st.session_state.get("pmv_expression") or {}
                st.caption(f"Loaded genes: {len(expr) if isinstance(expr, dict) else 0}")

                vmin, vmax = st.slider("Color scale range (log2FC)", min_value=-4.0, max_value=4.0, value=(float(st.session_state["pmv_vmin"]), float(st.session_state["pmv_vmax"])), step=0.1)
                st.session_state["pmv_vmin"] = float(vmin)
                st.session_state["pmv_vmax"] = float(vmax)
                cmap = st.selectbox("Colormap", options=["RdBu_r", "coolwarm", "viridis"], index=0, key="pmv_cmap")
                st.session_state["pmv_colormap"] = cmap

                if st.button("Apply Overlay", type="primary", key="pmv_apply_overlay"):
                    if not isinstance(expr, dict) or not expr:
                        st.warning("No expression data loaded.")
                    else:
                        try:
                            out = _api_post(
                                "/api/pathway-maps/overlay",
                                {"pathway_id": pid, "expression_data": expr, "colormap": cmap, "vmin": float(vmin), "vmax": float(vmax)},
                                timeout=120,
                            )
                            overlays = out.get("overlays") if isinstance(out, dict) else []
                            st.session_state["pmv_overlays"] = overlays if isinstance(overlays, list) else []
                        except Exception as e:  # noqa: BLE001
                            st.error(f"Overlay failed: {e}")

                render_color_legend(float(vmin), float(vmax), str(cmap))

    with tab4:
        st.subheader("Multi-Pathway")
        st.info("Multi-pathway comparison coming in Phase 2.")


__all__ = ["render_pathway_map_viewer_page"]


