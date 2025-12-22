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
SUM_ENDPOINT = f"{API_BASE}/api/notebook/notebook/summarize"


def _post_generate(query: str) -> Dict[str, Any]:
    with httpx.Client(timeout=120) as client:
        resp = client.post(GEN_ENDPOINT, json={"query": query})
        resp.raise_for_status()
        return resp.json()


def _post_summarize(notebook_json: Dict[str, Any]) -> Dict[str, Any]:
    with httpx.Client(timeout=120) as client:
        resp = client.post(SUM_ENDPOINT, json={"notebook_json": notebook_json})
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
            st.session_state.pop("generated_notebook_summary", None)
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

    col_a, col_b = st.columns([1, 2])
    with col_a:
        do_sum = st.button("Summarize", use_container_width=True)
    with col_b:
        st.caption("Uses the backend LLM to summarize this notebook.")

    if do_sum:
        try:
            with st.spinner("Summarizing notebook..."):
                st.session_state["generated_notebook_summary"] = _post_summarize(nb_json)
        except httpx.HTTPError as e:  # noqa: BLE001
            st.error(f"Failed to summarize notebook: {e}")

    summary = st.session_state.get("generated_notebook_summary")
    if isinstance(summary, dict) and summary:
        st.subheader("Summary")
        title = summary.get("title") or ""
        entity_summary = summary.get("entity_summary") or ""
        methods = summary.get("methods") or ""
        key_findings = summary.get("key_findings") or ""

        if title:
            st.markdown(f"**{title}**")
        if entity_summary:
            st.markdown("**Entity summary**")
            st.write(entity_summary)
        if methods:
            st.markdown("**Methods**")
            st.write(methods)
        if key_findings:
            st.markdown("**Key findings**")
            st.write(key_findings)

    st.subheader("Download")
    st.caption(f"Filename: {filename}")
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


