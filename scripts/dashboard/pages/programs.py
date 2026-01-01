"""Programs page for the Streamlit dashboard."""

from __future__ import annotations

import requests
import streamlit as st

from amprenta_rag.database.models import Program
from scripts.dashboard.db_session import db_session
from scripts.dashboard.pages.program_detail import render_pinned_dashboards_section
from scripts.dashboard.utils.cache import fetch_programs, clear_all_caches
from scripts.dashboard.utils.accessibility import (
    render_skip_link,
    add_heading_structure,
    ensure_minimum_contrast,
    accessible_text_input
)


def render_programs_page() -> None:
    """
    Render the Programs page with search and detailed views.

    Features:
    - Search programs by name (cached data)
    - Display program details
    - Show related datasets count
    - Cache refresh button
    """
    # Add accessibility features
    render_skip_link("main-programs-content")
    ensure_minimum_contrast()
    
    # Add main content landmark and heading
    st.markdown(
        """
        <main id="main-programs-content" role="main" aria-label="Research programs management">
        </main>
        """,
        unsafe_allow_html=True
    )
    
    add_heading_structure("🔬 Research Programs", level=1, id="programs-title")
    
    # Add cache refresh button
    col1, col2 = st.columns([10, 1])
    with col1:
        st.markdown("")  # Spacing since we added our own heading
    with col2:
        if st.button("🔄", help="Refresh data", key="refresh_programs"):
            clear_all_caches()
            st.rerun()

    # Use cached data
    programs_data = fetch_programs()
    
    search_term = st.text_input("Search programs by name", "")

    # Filter programs by search term (client-side filtering on cached data)
    if search_term:
        filtered_programs = [
            p for p in programs_data 
            if search_term.lower() in p.get("name", "").lower()
        ]
    else:
        filtered_programs = programs_data

    st.metric("Total Programs", len(filtered_programs))

    if filtered_programs:
            for program_data in filtered_programs:
                program_id = program_data.get("id")
                program_name = program_data.get("name", "Unknown")
                
                with st.expander(f"**{program_name}**"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{program_id}`")
                        if program_data.get("description"):
                            st.write(f"**Description:** {program_data['description']}")
                    with col2:
                        if program_data.get("created_at"):
                            st.write(f"**Created:** {program_data['created_at']}")
                        if program_data.get("disease"):
                            diseases = program_data["disease"]
                            if isinstance(diseases, list):
                                st.write(f"**Disease:** {', '.join(diseases)}")
                            else:
                                st.write(f"**Disease:** {diseases}")

                    # Show dataset count if available in cached data
                    if "dataset_count" in program_data:
                        st.write(f"**Related Datasets:** {program_data['dataset_count']}")

                    st.divider()
                    render_pinned_dashboards_section(program_id=program_id, key_prefix=f"prog_{program_id}_")

                    api_base = st.secrets.get("api_base_url", "http://localhost:8000")
                    if st.button("📄 Generate Report", key=f"report_{program_id}"):
                        payload = {
                            "entity_type": "program",
                            "entity_id": str(program_id),
                            "format": "html",
                        }
                        try:
                            resp = requests.post(
                                f"{api_base}/api/v1/reports/generate",
                                json=payload,
                                timeout=60,
                            )
                            if resp.ok:
                                data = resp.json()
                                link = data.get("download_url") or data.get("file_path")
                                if link:
                                    st.success("Report generated.")
                                    if isinstance(link, str) and link.startswith("http"):
                                        st.markdown(f"[Download report]({link})")
                                    else:
                                        st.code(link, language="text")
                                else:
                                    st.warning("Report generated but no link returned.")
                            else:
                                st.error(
                                    f"Failed to generate report ({resp.status_code}): {resp.text}"
                                )
                        except Exception as exc:  # pragma: no cover - UI fallback
                            st.error(f"Error generating report: {exc}")
    else:
        st.info("No programs found.")
