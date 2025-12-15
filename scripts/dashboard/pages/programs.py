"""Programs page for the Streamlit dashboard."""

from __future__ import annotations

import requests
import streamlit as st

from amprenta_rag.database.models import Program
from scripts.dashboard.db_session import db_session


def render_programs_page() -> None:
    """
    Render the Programs page with search and detailed views.

    Features:
    - Search programs by name
    - Display program details
    - Show related datasets count
    """
    st.header("🏥 Programs")

    with db_session() as db:
        search_term = st.text_input("Search programs by name", "")

        query = db.query(Program)
        if search_term:
            query = query.filter(Program.name.ilike(f"%{search_term}%"))

        programs = query.order_by(Program.created_at.desc()).all()

        st.metric("Total Programs", len(programs))

        if programs:
            for program in programs:
                with st.expander(f"**{program.name}**"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{program.id}`")
                        if program.description:
                            st.write(f"**Description:** {program.description}")
                    with col2:
                        st.write(f"**Created:** {program.created_at.strftime('%Y-%m-%d %H:%M')}")
                        if program.disease:
                            st.write(f"**Disease:** {', '.join(program.disease)}")

                    # Count related datasets (access while session is open)
                    dataset_count = len(program.datasets)
                    if dataset_count > 0:
                        st.write(f"**Related Datasets:** {dataset_count}")

                    api_base = st.secrets.get("api_base_url", "http://localhost:8000")
                    if st.button("📄 Generate Report", key=f"report_{program.id}"):
                        payload = {
                            "entity_type": "program",
                            "entity_id": str(program.id),
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
