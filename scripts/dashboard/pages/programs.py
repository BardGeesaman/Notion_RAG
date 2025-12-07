"""Programs page for the Streamlit dashboard."""

from __future__ import annotations

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
        else:
            st.info("No programs found.")
