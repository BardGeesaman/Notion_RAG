"""Connectivity Map dashboard page (LINCS / CMap)."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _list_signatures() -> List[Dict[str, Any]]:
    # Prefer v1 signatures API if available.
    with httpx.Client(timeout=30) as client:
        r = client.get(f"{API_BASE}/api/v1/signatures")
        r.raise_for_status()
        out = r.json()
        return list(out or [])


def render_connectivity_map_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("🧬 Connectivity Map")
    st.caption("Compute CMap connectivity between internal signatures and LINCS (L1000) signatures.")

    # Load internal signatures (best-effort; page still renders if API fails)
    signatures: List[Dict[str, Any]] = []
    sig_error: Optional[str] = None
    try:
        signatures = _list_signatures()
    except Exception as e:  # noqa: BLE001
        sig_error = str(e)

    sig_map: Dict[str, str] = {}
    for s in signatures:
        sid = s.get("id") or s.get("signature_id")  # some endpoints use "id"
        name = s.get("name") or s.get("short_id") or str(sid)
        if sid:
            sig_map[f"{name} ({str(sid)[:8]})"] = str(sid)

    tab1, tab2, tab3 = st.tabs(["Signature", "Top Reversals", "Top Mimics"])

    with tab1:
        st.subheader("Signature Selector")
        if sig_error:
            st.error(f"Failed to load signatures: {sig_error}")
        if not sig_map:
            st.selectbox("Internal Signature", ["(no signatures available)"])
            st.button("Compute Connectivity", disabled=True)
        else:
            sel = st.selectbox("Internal Signature", list(sig_map.keys()))
            signature_id = sig_map[sel]
            if st.button("Compute Connectivity", type="primary"):
                try:
                    with httpx.Client(timeout=300) as client:
                        rr = client.post(f"{API_BASE}/api/connectivity/compute/{signature_id}")
                    if rr.status_code >= 400:
                        st.error(f"Compute failed ({rr.status_code}): {rr.text}")
                    else:
                        st.success("Connectivity computed.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Request failed: {e}")

            st.session_state["connectivity_selected_signature_id"] = signature_id

    selected_id = st.session_state.get("connectivity_selected_signature_id")

    with tab2:
        st.subheader("Top Reversals")
        threshold = st.slider("Score threshold (reversals)", min_value=-1.0, max_value=0.0, value=-0.1, step=0.01)
        if not selected_id:
            st.info("Select a signature in the Signature tab first.")
        else:
            try:
                with httpx.Client(timeout=30) as client:
                    r = client.get(f"{API_BASE}/api/connectivity/signatures/{selected_id}/top-reversals", params={"n": 200})
                r.raise_for_status()
                rows = [x for x in (r.json() or []) if float(x.get("score", 0.0)) <= float(threshold)]
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load reversals: {e}")
                rows = []

            if not rows:
                st.info("No reversals found (or scores not computed yet).")
            else:
                df = pd.DataFrame(
                    [
                        {
                            "rank": i + 1,
                            "compound": r.get("pert_iname"),
                            "pert_id": r.get("pert_id"),
                            "cell": r.get("cell_id"),
                            "score": r.get("score"),
                        }
                        for i, r in enumerate(rows)
                    ]
                )
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab3:
        st.subheader("Top Mimics")
        threshold = st.slider("Score threshold (mimics)", min_value=0.0, max_value=1.0, value=0.1, step=0.01)
        if not selected_id:
            st.info("Select a signature in the Signature tab first.")
        else:
            try:
                with httpx.Client(timeout=30) as client:
                    r = client.get(f"{API_BASE}/api/connectivity/signatures/{selected_id}/top-mimics", params={"n": 200})
                r.raise_for_status()
                rows = [x for x in (r.json() or []) if float(x.get("score", 0.0)) >= float(threshold)]
            except Exception as e:  # noqa: BLE001
                st.error(f"Failed to load mimics: {e}")
                rows = []

            if not rows:
                st.info("No mimics found (or scores not computed yet).")
            else:
                df = pd.DataFrame(
                    [
                        {
                            "rank": i + 1,
                            "compound": r.get("pert_iname"),
                            "pert_id": r.get("pert_id"),
                            "cell": r.get("cell_id"),
                            "score": r.get("score"),
                        }
                        for i, r in enumerate(rows)
                    ]
                )
                st.dataframe(df, use_container_width=True, hide_index=True)


__all__ = ["render_connectivity_map_page"]


