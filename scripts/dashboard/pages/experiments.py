"""Experiments page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Experiment
from scripts.dashboard.db_session import db_session


def render_experiments_page() -> None:
    """
    Render the Experiments page with search, summary, and detailed views.

    Features:
    - Search experiments by name
    - Summary table with export
    - Detailed experiment views
    - Show related datasets count
    """
    st.header("🔬 Experiments")

    with db_session() as db:
        search_term = st.text_input("Search experiments by name", "")

        query = db.query(Experiment)
        if search_term:
            query = query.filter(Experiment.name.ilike(f"%{search_term}%"))

        experiments = query.order_by(Experiment.created_at.desc()).all()

        st.metric("Total Experiments", len(experiments))

        if experiments:
            # Summary table
            experiment_data = []
            for exp in experiments:
                experiment_data.append(
                    {
                        "Name": exp.name,
                        "Type": exp.type or "",
                        "Disease": ", ".join(exp.disease) if exp.disease else "",
                        "Matrix": ", ".join(exp.matrix) if exp.matrix else "",
                        "Created": exp.created_at.strftime("%Y-%m-%d"),
                        "Datasets": len(exp.datasets),
                    }
                )
            df_experiments = pd.DataFrame(experiment_data)
            st.dataframe(df_experiments, use_container_width=True, hide_index=True)

            # Export button
            csv_experiments = df_experiments.to_csv(index=False)
            st.download_button(
                label="📥 Download Experiments (CSV)",
                data=csv_experiments,
                file_name="experiments.csv",
                mime="text/csv",
            )

            st.markdown("---")
            st.subheader("Experiment Details")

            for experiment in experiments:
                with st.expander(f"**{experiment.name}**"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{experiment.id}`")
                        if experiment.type:
                            st.write(f"**Type:** {experiment.type}")
                        if experiment.description:
                            st.write(f"**Description:** {experiment.description}")
                    with col2:
                        st.write(f"**Created:** {experiment.created_at.strftime('%Y-%m-%d %H:%M')}")
                        if experiment.disease:
                            st.write(f"**Disease:** {', '.join(experiment.disease)}")
                        if experiment.matrix:
                            st.write(f"**Matrix:** {', '.join(experiment.matrix)}")

                    # Count related datasets
                    dataset_count = len(experiment.datasets)
                    if dataset_count > 0:
                        st.write(f"**Related Datasets:** {dataset_count}")
        else:
            st.info("No experiments found.")
