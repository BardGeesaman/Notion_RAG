"""Experiments page for the Streamlit dashboard."""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Experiment
from scripts.dashboard.db_session import db_session


def _render_browse_tab() -> None:
    with db_session() as db:
        search_term = st.text_input("Search experiments by name", "", key="exp_search")

        query = db.query(Experiment)
        if search_term:
            query = query.filter(Experiment.name.ilike(f"%{search_term}%"))

        experiments: List[Experiment] = query.order_by(Experiment.created_at.desc()).all()

        st.metric("Total Experiments", len(experiments))

        if experiments:
            # Summary table
            experiment_data: List[Dict[str, Any]] = []
            for exp in experiments:
                experiment_data.append(
                    {
                        "Name": exp.name,
                        "Type": exp.type or "",
                        "Design": exp.design_type or "",
                        "Disease": ", ".join(exp.disease) if exp.disease else "",
                        "Matrix": ", ".join(exp.matrix) if exp.matrix else "",
                        "Created": exp.created_at.strftime("%Y-%m-%d"),
                        "Datasets": len(exp.datasets),
                    }
                )
            df_experiments = pd.DataFrame(experiment_data)
            st.dataframe(df_experiments, width='stretch', hide_index=True)

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
                        st.write(f"**Design:** {experiment.design_type or 'n/a'}")
                    with col2:
                        st.write(f"**Created:** {experiment.created_at.strftime('%Y-%m-%d %H:%M')}")
                        if experiment.disease:
                            st.write(f"**Disease:** {', '.join(experiment.disease)}")
                        if experiment.matrix:
                            st.write(f"**Matrix:** {', '.join(experiment.matrix)}")
                        if experiment.sample_groups:
                            st.write("**Sample Groups:**")
                            st.json(experiment.sample_groups)

                    dataset_count = len(experiment.datasets)
                    if dataset_count > 0:
                        st.write(f"**Related Datasets:** {dataset_count}")
        else:
            st.info("No experiments found.")


def _render_edit_tab() -> None:
    st.subheader("Edit Experimental Design")
    with db_session() as db:
        experiments: List[Experiment] = db.query(Experiment).order_by(Experiment.created_at.desc()).all()

        st.markdown("### Batch Operations")
        if st.button("Auto-detect Design Types", type="secondary"):
            from amprenta_rag.ingestion.design_integration import batch_apply_design_extraction

            with st.spinner("Detecting design types..."):
                res = batch_apply_design_extraction()
            if "error" in res:
                st.error(f"Auto-detection failed: {res['error']}")
            else:
                st.success(
                    f"Auto-detection complete: {res.get('updated', 0)} updated out of {res.get('processed', 0)} processed."
                )
        st.markdown("---")
        if not experiments:
            st.info("No experiments available to edit.")
            return

        exp_options = {exp.name: exp for exp in experiments}
        selected_name = st.selectbox("Select experiment", list(exp_options.keys()))
        experiment = exp_options[selected_name]

        design_type = st.selectbox(
            "Design type",
            [
                "case_control",
                "time_course",
                "intervention",
                "dose_response",
                "multi_factorial",
                "observational",
            ],
            index=0 if not experiment.design_type else  \
                [
                    "case_control",
                    "time_course",
                    "intervention",
                    "dose_response",
                    "multi_factorial",
                    "observational",
                ].index(experiment.design_type)
                if experiment.design_type in [
                    "case_control",
                    "time_course",
                    "intervention",
                    "dose_response",
                    "multi_factorial",
                    "observational",
                ]
                else 0,
        )

        sample_groups_text = st.text_area(
            "Sample groups (JSON)",
            value=json.dumps(experiment.sample_groups or {}, indent=2),
            height=160,
        )

        timepoints_text = st.text_input(
            "Timepoints (comma-separated, optional)",
            value=", ".join(experiment.design_metadata.get("timepoints", [])) if experiment.design_metadata else "",
        )

        if st.button("Save design", type="primary"):
            try:
                sample_groups = json.loads(sample_groups_text) if sample_groups_text.strip() else {}
            except Exception as exc:
                st.error(f"Invalid sample groups JSON: {exc}")
                return

            metadata: Dict[str, Any] = experiment.design_metadata or {}
            if timepoints_text.strip():
                metadata["timepoints"] = [tp.strip() for tp in timepoints_text.split(",") if tp.strip()]

            experiment.design_type = design_type
            experiment.sample_groups = sample_groups
            experiment.design_metadata = metadata

            db.add(experiment)
            db.commit()
            st.success("Design updated.")


def render_experiments_page() -> None:
    """
    Render the Experiments page with search, summary, and design editing.
    """
    st.header("🔬 Experiments")

    tabs = st.tabs(["Browse Experiments", "Edit Design"])
    with tabs[0]:
        _render_browse_tab()
    with tabs[1]:
        _render_edit_tab()


__all__ = ["render_experiments_page"]
