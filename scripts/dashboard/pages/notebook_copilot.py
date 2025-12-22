"""Notebook Co-Pilot page for the Streamlit dashboard."""

from __future__ import annotations

import streamlit as st

from scripts.dashboard.components.copilot_sidebar import render_copilot_sidebar


def render_notebook_copilot_page() -> None:
    """Render the Notebook Co-Pilot page."""
    st.header("🧠 Notebook Co-Pilot")
    st.caption("Generate notebook-ready Python cells from context + intent.")

    with st.expander("How to use", expanded=True):
        st.markdown(
            """
- **Pick an action** (Load Dataset / HTS QC / Dose-response / Publish)
- **Select the entity** (dataset/campaign/compound/etc.)
- Optionally add a **free-form prompt** for more detail
- Click **Generate code** and copy/paste into your notebook

Notes:
- Copilot requires LLM credentials (e.g. OpenAI) in the environment.
- You can set the model via `AMPRENTA_NOTEBOOK_MODEL` (default: `gpt-4o-mini`).
            """.strip()
        )

    st.markdown("---")

    # Render copilot UI in the main page body (not the sidebar).
    render_copilot_sidebar()


