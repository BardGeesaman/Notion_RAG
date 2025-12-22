"""Experiments page for the Streamlit dashboard."""

from __future__ import annotations

import json
from typing import Any, Dict, List
from uuid import UUID

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Experiment, ExperimentTemplate
from amprenta_rag.auth.session import get_current_user
from amprenta_rag.ingestion.design_integration import (
    batch_apply_design_extraction,
    create_experiment_from_study,
)
from amprenta_rag.ingestion.repositories import fetch_study_metadata
from amprenta_rag.utils.data_export import export_experiments
from amprenta_rag.utils.saved_filters import save_filter, get_user_filters, delete_filter
from amprenta_rag.analysis.study_critique import assess_study_quality
from amprenta_rag.notifications.email_service import send_share_notification, is_email_configured
from amprenta_rag.export.slide_generator import generate_experiment_slides
from amprenta_rag.utils.optimistic_lock import update_with_lock, ConflictError
from scripts.dashboard.db_session import db_session


def _handle_reload(exp_id: int, actual_version: int) -> None:
    """Callback for Reload Latest button - updates version in session state."""
    st.session_state[f"edit_version_{exp_id}"] = actual_version
    st.session_state.pop(f"conflict_{exp_id}", None)


def _render_browse_tab() -> None:
    with db_session() as db:
        # Saved Filters section
        user = get_current_user()
        if user:
            saved_filters = get_user_filters("experiment", user.get("id"), db)

            col1, col2, col3 = st.columns([2, 1, 1])
            with col1:
                if saved_filters:
                    filter_names = [f.name for f in saved_filters]
                    selected_filter_name = st.selectbox(
                        "📌 Saved Filters",
                        ["None"] + filter_names,
                        key="saved_filter_select"
                    )
                    if selected_filter_name != "None":
                        selected_filter = next(f for f in saved_filters if f.name == selected_filter_name)
                        st.info(f"Filter '{selected_filter_name}' selected. Filters: {selected_filter.filters}")
            with col2:
                if st.session_state.get("show_save_filter", False):
                    with st.form("save_filter_form"):
                        filter_name = st.text_input("Filter name", key="filter_name_input")
                        submitted = st.form_submit_button("Save", type="primary")
                        if submitted and filter_name:
                            # Get current filter state (would need to capture from UI)
                            current_filters = {"search": st.session_state.get("exp_search", "")}  # Capture search term
                            save_filter(filter_name, "experiment", current_filters, user.get("id"), db)
                            st.session_state["show_save_filter"] = False
                            st.success(f"Filter '{filter_name}' saved!")
                            st.rerun()
                        if st.form_submit_button("Cancel"):
                            st.session_state["show_save_filter"] = False
                            st.rerun()
                else:
                    if st.button("💾 Save Current", key="save_filter_btn"):
                        st.session_state["show_save_filter"] = True
                        st.rerun()
            with col3:
                if saved_filters and st.session_state.get("saved_filter_select") != "None":
                    selected_filter_name = st.session_state.get("saved_filter_select")
                    if selected_filter_name and selected_filter_name != "None":
                        if st.button("🗑️ Delete", key="delete_filter_btn"):
                            selected_filter = next(f for f in saved_filters if f.name == selected_filter_name)
                            delete_filter(str(selected_filter.id), db)
                            st.success("Filter deleted!")
                            st.rerun()

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

            # Export section
            st.markdown("### Export")
            col1, col2 = st.columns([1, 3])
            with col1:
                export_format = st.selectbox("Format", ["csv", "json", "excel"], key="exp_export_format")
            with col2:
                experiment_ids = [str(exp.id) for exp in experiments]
                mime_types = {"csv": "text/csv", "json": "application/json", "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}
                file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                try:
                    export_data = export_experiments(experiment_ids, export_format, db)
                    st.download_button(
                        label=f"📥 Download Experiments ({export_format.upper()})",
                        data=export_data,
                        file_name=f"experiments.{file_extensions[export_format]}",
                        mime=mime_types[export_format],
                    )
                except Exception as e:
                    st.error(f"Export failed: {e}")

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

                    # Export as PowerPoint button
                    try:
                        pptx_data = generate_experiment_slides(experiment.id, db)
                        experiment_name_safe = (experiment.name or "experiment").replace(" ", "_").replace("/", "_")[:50]
                        st.download_button(
                            label="📊 Export as PowerPoint",
                            data=pptx_data,
                            file_name=f"{experiment_name_safe}.pptx",
                            mime_type="application/vnd.openxmlformats-officedocument.presentationml.presentation",
                            key=f"export_pptx_{experiment.id}",
                        )
                    except Exception as e:
                        st.error(f"Failed to generate PowerPoint: {e}")

                    # Share via Email section
                    if is_email_configured():
                        st.markdown("---")
                        st.markdown("### 📧 Share via Email")
                        with st.form(f"share_email_form_{experiment.id}"):
                            recipient_email = st.text_input("Recipient Email", key=f"recipient_email_{experiment.id}")
                            share_message = st.text_area("Optional Message", key=f"share_message_{experiment.id}", placeholder="Add a personal message...")
                            if st.form_submit_button("Send", type="primary"):
                                if recipient_email:
                                    user = get_current_user()
                                    from_user = user.get("username", "Unknown User") if user else "Unknown User"
                                    try:
                                        success = send_share_notification(
                                            entity_type="Experiment",
                                            entity_id=experiment.id,
                                            to=recipient_email,
                                            from_user=from_user,
                                            message=share_message if share_message.strip() else None,
                                            db=db,
                                        )
                                        if success:
                                            st.success(f"Email sent successfully to {recipient_email}!")
                                        else:
                                            st.error("Failed to send email. Please check the logs.")
                                    except Exception as e:
                                        st.error(f"Error sending email: {e}")
                                else:
                                    st.error("Please enter a recipient email address")

                    # Save as Template button
                    user = get_current_user()
                    if user:
                        if st.button("💾 Save as Template", key=f"save_template_{experiment.id}"):
                            with st.form(f"template_form_{experiment.id}"):
                                template_name = st.text_input("Template Name*", value=f"{experiment.name} Template")
                                template_description = st.text_area("Description", value=experiment.description or "")

                                if st.form_submit_button("Create Template", type="primary"):
                                    if template_name:
                                        try:
                                            template = ExperimentTemplate(
                                                name=template_name,
                                                description=template_description if template_description.strip() else None,
                                                design_type=experiment.design_type,
                                                organism=None,  # Can be added later if needed
                                                sample_groups=experiment.sample_groups,
                                                created_by_id=UUID(user.get("id")) if user.get("id") and user.get("id") != "test" else None,
                                            )
                                            db.add(template)
                                            db.commit()
                                            st.success(f"Template '{template_name}' created!")
                                            st.rerun()
                                        except Exception as e:
                                            st.error(f"Failed to create template: {e}")
                                    else:
                                        st.error("Template name is required")
        else:
            st.info("No experiments found.")


def _render_edit_tab() -> None:
    st.subheader("Edit Experimental Design")
    with db_session() as db:
        experiments: List[Experiment] = db.query(Experiment).order_by(Experiment.created_at.desc()).all()

        st.markdown("### Batch Operations")
        if st.button("Auto-detect Design Types", type="secondary"):
            with st.spinner("Detecting design types..."):
                res = batch_apply_design_extraction()
            if "error" in res:
                st.error(f"Auto-detection failed: {res['error']}")
            else:
                st.success(
                    f"Auto-detection complete: {res.get('updated', 0)} updated out of {res.get('processed', 0)} processed."
                )
        st.markdown("---")

        st.markdown("### Import from Repository")
        repo_options = ["GEO", "ArrayExpress", "MW", "MW_METABOLOMICS", "MW_LIPIDOMICS", "PRIDE", "MetaboLights", "ENA"]
        repo = st.selectbox("Repository", repo_options, index=0)
        study_id = st.text_input("Study ID (e.g., GSE12345)")
        if st.button("Create Experiment from Study"):
            with st.spinner("Fetching study and creating experiment..."):
                try:
                    md = fetch_study_metadata(study_id, repository=repo)
                except Exception as exc:
                    md = None
                    st.error(f"Failed to fetch study metadata: {exc}")
                if md:
                    exp_id = create_experiment_from_study(md)
                    if exp_id:
                        st.success(f"Experiment created: {exp_id}")
                    else:
                        st.error("Failed to create experiment from study.")
                elif md is not None:
                    st.error("No metadata returned for that study.")
        st.markdown("---")

        if not experiments:
            st.info("No experiments available to edit.")
            return

        exp_options = {exp.name: exp for exp in experiments}
        selected_name = st.selectbox("Select experiment", list(exp_options.keys()))
        experiment = exp_options[selected_name]

        # Store version for optimistic locking (only on first load, not reruns)
        if f"edit_version_{experiment.id}" not in st.session_state:
            st.session_state[f"edit_version_{experiment.id}"] = experiment.version

        # Debug: show current version tracking
        st.caption(
            f"[DEBUG] exp_id={experiment.id}, edit_version={st.session_state.get(f'edit_version_{experiment.id}', 'NOT SET')}, "
            f"db_version={experiment.version}"
        )
        conflict_key = f"conflict_{experiment.id}"
        if conflict_key in st.session_state:
            st.caption(f"[DEBUG] conflict_state={st.session_state[conflict_key]}")

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
            index=0
            if not experiment.design_type
            else [
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

            updates = {
                "design_type": design_type,
                "sample_groups": sample_groups,
                "design_metadata": metadata,
            }

            expected_version = st.session_state.get(f"edit_version_{experiment.id}", 1)

            try:
                update_with_lock(experiment, updates, expected_version, db)
                st.success("Saved successfully!")
                # Update stored version
                st.session_state[f"edit_version_{experiment.id}"] = experiment.version
            except ConflictError as e:
                st.session_state[f"conflict_{experiment.id}"] = {"expected": e.expected, "actual": e.actual}
                st.rerun()

        conflict_key = f"conflict_{experiment.id}"
        if conflict_key in st.session_state:
            conflict = st.session_state[conflict_key]
            st.error("⚠️ Conflict: This record was modified by another user.")
            st.warning(f"Your version: {conflict['expected']}, Current version: {conflict['actual']}")
            col1, col2 = st.columns(2)
            with col1:
                st.button(
                    "Reload Latest",
                    key=f"reload_{experiment.id}",
                    on_click=_handle_reload,
                    args=(experiment.id, conflict["actual"]),
                )
            with col2:
                st.button(
                    "Force Save",
                    key=f"force_{experiment.id}",
                    on_click=_handle_reload,
                    args=(experiment.id, conflict["actual"]),
                )


def _render_templates_tab() -> None:
    """Render the Templates tab."""
    user = get_current_user()
    user_id = user.get("id") if user else None

    with db_session() as db:
        templates = db.query(ExperimentTemplate).order_by(ExperimentTemplate.created_at.desc()).all()

        st.subheader("Experiment Templates")
        st.markdown("Use templates to quickly create experiments with pre-filled values.")

        if templates:
            for template in templates:
                with st.expander(f"**{template.name}**"):
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.write(template.description or "_No description_")
                        if template.design_type:
                            st.write(f"**Design Type:** {template.design_type}")
                        if template.organism:
                            st.write(f"**Organism:** {', '.join(template.organism)}")
                        if template.sample_groups:
                            st.write("**Sample Groups:**")
                            st.json(template.sample_groups)
                    with col2:
                        if st.button("Use Template", key=f"use_template_{template.id}"):
                            st.session_state["use_template_id"] = str(template.id)
                            st.session_state["template_data"] = {
                                "name": template.name,
                                "description": template.description,
                                "design_type": template.design_type,
                                "organism": template.organism,
                                "sample_groups": template.sample_groups,
                            }
                            st.rerun()
        else:
            st.info("No templates available. Create one by saving an experiment as a template.")

        st.markdown("---")
        st.subheader("Create from Template")

        if st.session_state.get("template_data"):
            template_data = st.session_state["template_data"]
            st.info("Template loaded. Fill in the form below and create your experiment.")

            with st.form("create_from_template"):
                name = st.text_input("Experiment Name*", value=template_data.get("name", ""))
                description = st.text_area("Description", value=template_data.get("description", ""))
                design_type = st.selectbox(
                    "Design Type",
                    ["case_control", "time_course", "intervention", "dose_response", "multi_factorial", "observational"],
                    index=["case_control", "time_course", "intervention", "dose_response", "multi_factorial", "observational"].index(template_data.get("design_type", "case_control")) if template_data.get("design_type") in ["case_control", "time_course", "intervention", "dose_response", "multi_factorial", "observational"] else 0,
                )
                organism_input = st.text_input("Organism (comma-separated)", value=", ".join(template_data.get("organism", [])) if template_data.get("organism") else "")
                sample_groups_json = st.text_area("Sample Groups (JSON)", value=json.dumps(template_data.get("sample_groups", {}), indent=2), height=100)

                if st.form_submit_button("Create Experiment", type="primary"):
                    if not name:
                        st.error("Experiment name is required")
                    else:
                        try:
                            sample_groups = json.loads(sample_groups_json) if sample_groups_json.strip() else {}
                            organism = [o.strip() for o in organism_input.split(",") if o.strip()] if organism_input.strip() else None

                            experiment = Experiment(
                                name=name,
                                description=description if description.strip() else None,
                                design_type=design_type,
                                organism=organism,
                                sample_groups=sample_groups,
                                created_by_id=UUID(user_id) if user_id and user_id != "test" else None,
                            )

                            db.add(experiment)
                            db.commit()
                            db.refresh(experiment)

                            # Fire workflow trigger
                            from amprenta_rag.automation.engine import fire_trigger
                            fire_trigger("experiment_created", {
                                "experiment_id": str(experiment.id),
                                "name": experiment.name
                            }, db)

                            st.success(f"Experiment '{name}' created successfully!")
                            st.session_state.pop("template_data", None)
                            st.session_state.pop("use_template_id", None)
                            st.rerun()
                        except json.JSONDecodeError as e:
                            st.error(f"Invalid JSON in sample groups: {e}")
                        except Exception as e:
                            st.error(f"Failed to create experiment: {e}")


def render_experiments_page() -> None:
    """Render the Experiments page with search, summary, and design editing."""
    st.header("🔬 Experiments")

    tabs = st.tabs(["Browse Experiments", "Edit Design", "Templates", "Study Quality"])
    with tabs[0]:
        _render_browse_tab()
    with tabs[1]:
        _render_edit_tab()
    with tabs[2]:
        _render_templates_tab()
    with tabs[3]:
        _render_study_quality_tab()


def _render_study_quality_tab() -> None:
    """Render the Study Quality tab."""
    st.subheader("Study Quality Assessment")
    st.markdown("Automated quality assessment for experiments.")

    with db_session() as db:
        experiments = db.query(Experiment).order_by(Experiment.created_at.desc()).limit(100).all()

        if not experiments:
            st.info("No experiments available for assessment.")
            return

        exp_options = {exp.name: exp.id for exp in experiments}
        selected_name = st.selectbox("Select Experiment", list(exp_options.keys()))

        if st.button("🔍 Run Quality Assessment", type="primary"):
            exp_id = exp_options[selected_name]
            with st.spinner("Assessing study quality..."):
                result = assess_study_quality(exp_id, db)
                st.session_state["quality_assessment"] = result
                st.session_state["assessed_experiment_id"] = str(exp_id)
                st.rerun()

        if st.session_state.get("quality_assessment") and st.session_state.get("assessed_experiment_id") == str(exp_options.get(selected_name)):
            result = st.session_state["quality_assessment"]
            quality_score = result.get("quality_score", 0)

            # Quality Score with color
            if quality_score >= 80:
                delta_color = "normal"
                st.metric("Quality Score", f"{quality_score}/100", delta=None, delta_color=delta_color)
                st.success("✅ High Quality")
            elif quality_score >= 50:
                delta_color = "off"
                st.metric("Quality Score", f"{quality_score}/100", delta=None, delta_color=delta_color)
                st.warning("⚠️ Moderate Quality")
            else:
                delta_color = "inverse"
                st.metric("Quality Score", f"{quality_score}/100", delta=None, delta_color=delta_color)
                st.error("❌ Low Quality")

            st.markdown("---")

            # Overall Summary
            st.markdown("### Summary")
            st.info(result.get("summary", "No summary available"))

            st.markdown("---")

            # Design Flaws
            issues = result.get("issues", [])
            design_flaws = [i for i in issues if "flaw" in i]
            data_gaps = [i for i in issues if i.get("flaw") in ["Missing organism field", "Missing description", "No linked datasets"]]

            if design_flaws:
                st.markdown("### ⚠️ Design Flaws")
                flaw_data = []
                for flaw in design_flaws:
                    severity_emoji = {
                        "high": "🔴",
                        "medium": "🟡",
                        "low": "🟢"
                    }.get(flaw.get("severity", "medium"), "⚪")

                    flaw_data.append({
                        "Severity": f"{severity_emoji} {flaw.get('severity', 'medium').upper()}",
                        "Issue": flaw.get("flaw", ""),
                        "Suggestion": flaw.get("suggestion", ""),
                    })

                df_flaws = pd.DataFrame(flaw_data)
                st.dataframe(df_flaws, use_container_width=True, hide_index=True)

            # Data Gaps
            if data_gaps:
                st.markdown("### 📋 Data Gaps")
                for gap in data_gaps:
                    st.markdown(f"- {gap.get('flaw', '')}")

            if not design_flaws and not data_gaps:
                st.success("✅ No issues detected! Study quality is good.")


__all__ = ["render_experiments_page"]
