"""Query→Notebook Generator page for the Streamlit dashboard."""

from __future__ import annotations

import json
import os
from typing import Any, Dict, Optional

import httpx
import nbformat
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")
GEN_ENDPOINT = f"{API_BASE}/api/notebook/notebook/generate-from-query"


def _post_generate(query: str) -> Dict[str, Any]:
    with httpx.Client(timeout=120) as client:
        resp = client.post(GEN_ENDPOINT, json={"query": query})
        resp.raise_for_status()
        return resp.json()


def render_notebook_generator_page() -> None:
    """Render the Query→Notebook Generator page."""
    st.header("🧾 Query → Notebook Generator")
    st.caption("Ask a question and generate a full Jupyter notebook (.ipynb).")

    query = st.text_area("Question", height=120, placeholder="e.g., Summarize datasets related to ALS ceramides and suggest follow-up analyses")

    col1, col2 = st.columns([1, 2])
    with col1:
        do_gen = st.button("Generate Notebook", type="primary", disabled=not query.strip())
    with col2:
        st.caption("Requires API + LLM credentials configured on the backend.")

    if do_gen:
        try:
            with st.spinner("Generating notebook..."):
                payload = _post_generate(query.strip())
            st.session_state["generated_notebook"] = payload
        except httpx.HTTPError as e:  # noqa: BLE001
            st.error(f"Failed to generate notebook: {e}")
            return

    payload = st.session_state.get("generated_notebook")
    if not payload:
        return

    nb_json = payload.get("notebook_json")
    filename = payload.get("filename") or "generated.ipynb"

    if not isinstance(nb_json, dict):
        st.error("Notebook payload invalid.")
        return

    st.subheader("Preview")
    nb = nbformat.from_dict(nb_json)
    for i, cell in enumerate(nb.get("cells", []), start=1):
        ctype = cell.get("cell_type")
        src = cell.get("source") or ""
        with st.expander(f"Cell {i} ({ctype})", expanded=(i <= 3)):
            if ctype == "markdown":
                st.markdown(src)
            else:
                st.code(src, language="python")

    st.subheader("Download")
    nb_bytes = nbformat.writes(nb).encode("utf-8")
    st.download_button(
        "Download .ipynb",
        data=nb_bytes,
        file_name=filename,
        mime="application/x-ipynb+json",
        use_container_width=True,
    )

    st.caption("Raw notebook JSON")
    st.code(json.dumps(nb_json, indent=2)[:4000], language="json")


__all__ = ["render_notebook_generator_page"]


