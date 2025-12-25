"""Notebook template gallery (Jupyter Advanced)."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import streamlit as st

from scripts.dashboard.components.notebook_card import render_notebook_card
from scripts.dashboard.components.review_card import render_review_badge


API_BASE = os.environ.get("API_URL", "http://localhost:8000")
JUPYTERHUB_URL = os.environ.get("JUPYTERHUB_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _score(query: str, text: str) -> float:
    q = (query or "").strip().lower()
    t = (text or "").strip().lower()
    if not q:
        return 1.0
    if q in t:
        return 1.0
    # simple fuzzy fallback: overlap ratio
    q_tokens = set(q.split())
    t_tokens = set(t.split())
    if not q_tokens:
        return 0.0
    return len(q_tokens & t_tokens) / len(q_tokens)


def _truncate(s: str, n: int = 160) -> str:
    t = (s or "").strip()
    if len(t) <= n:
        return t
    return t[: n - 1].rstrip() + "…"


@st.cache_data(ttl=30, show_spinner=False)
def _review_status(nb_path: str) -> dict:
    try:
        return dict(_api_get(f"/api/reviews/notebook/{nb_path}", timeout=10) or {})
    except Exception:
        return {}


def render_notebook_gallery_page() -> None:
    from scripts.dashboard.auth import require_auth

    require_auth()

    st.header("Notebook Gallery")
    st.caption("Browse curated Jupyter templates and open them in JupyterHub.")

    # Controls
    query = st.text_input("Search", placeholder="Search title/description/tags")

    try:
        templates: List[Dict[str, Any]] = list(_api_get("/api/notebooks/templates") or [])
    except Exception as e:  # noqa: BLE001
        st.error(f"Failed to load templates: {e}")
        return

    all_tags = sorted({t for tpl in templates for t in (tpl.get("tags") or []) if t})
    all_omics = sorted({str(tpl.get("omics_type") or "").strip() for tpl in templates if tpl.get("omics_type")})
    all_omics = [o for o in all_omics if o]

    c1, c2 = st.columns(2)
    with c1:
        tags = st.multiselect("Tags", options=all_tags, default=[])
    with c2:
        omics_type = st.selectbox("Omics type", options=["(any)"] + all_omics, index=0)

    # Filter
    tagset = set(tags)
    filtered: List[Dict[str, Any]] = []
    for tpl in templates:
        if omics_type != "(any)" and str(tpl.get("omics_type") or "").strip() != omics_type:
            continue
        if tagset:
            tset = set(tpl.get("tags") or [])
            if not tagset.issubset(tset):
                continue
        blob = " ".join(
            [
                str(tpl.get("title") or ""),
                str(tpl.get("description") or ""),
                " ".join([str(x) for x in (tpl.get("tags") or [])]),
            ]
        )
        if _score(query, blob) <= 0:
            continue
        filtered.append(tpl)

    # Sort by best score for query
    filtered.sort(
        key=lambda tpl: _score(
            query,
            " ".join([str(tpl.get("title") or ""), str(tpl.get("description") or ""), " ".join(tpl.get("tags") or [])]),
        ),
        reverse=True,
    )

    st.write(f"Showing {len(filtered)} template(s).")

    cols = st.columns(3)
    for i, tpl in enumerate(filtered):
        col = cols[i % 3]
        with col:
            tpl_id = str(tpl.get("id"))
            nb_path = str(tpl.get("notebook_path") or "")

            source_url = f"{API_BASE}/api/notebooks/templates/{tpl_id}/download" if tpl_id else None
            if nb_path:
                render_review_badge(_review_status(nb_path))
            render_notebook_card(
                notebook_path=nb_path,
                title=str(tpl.get("title") or tpl.get("id") or "Untitled"),
                description=_truncate(str(tpl.get("description") or "")),
                tags=list(tpl.get("tags") or []),
                # JupyterHub link preserved (Open in Jupyter)
                jupyter_url=f"{JUPYTERHUB_URL.rstrip('/')}/notebooks/{nb_path}" if nb_path else None,
                # Voila: /voila/render/{notebook_path}
                voila_url=f"{JUPYTERHUB_URL.rstrip('/')}/voila/render/{nb_path}" if nb_path else None,
                source_url=source_url,
                preview_image=str(tpl.get("preview_image") or "") if tpl.get("preview_image") else None,
            )

            # Download button (fetch from API; keep existing behavior)
            try:
                with httpx.Client(timeout=30) as client:
                    r = client.get(f"{API_BASE}/api/notebooks/templates/{tpl_id}/download")
                if r.status_code == 200:
                    st.download_button(
                        "Download .ipynb",
                        data=r.content,
                        file_name=nb_path.split("/")[-1] if nb_path else f"{tpl_id}.ipynb",
                        mime="application/x-ipynb+json",
                        key=f"dl_{tpl_id}",
                    )
                else:
                    st.button("Download .ipynb", disabled=True, key=f"dl_disabled_{tpl_id}")
            except Exception:
                st.button("Download .ipynb", disabled=True, key=f"dl_err_{tpl_id}")


__all__ = ["render_notebook_gallery_page"]


