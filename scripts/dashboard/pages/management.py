"""Data Management page for the Streamlit dashboard."""

from __future__ import annotations

import streamlit as st

from amprenta_rag.api.schemas import DatasetUpdate, ExperimentUpdate, ProgramUpdate
from amprenta_rag.api.services.datasets import update_dataset
from amprenta_rag.api.services.experiments import update_experiment
from amprenta_rag.api.services.programs import update_program
from amprenta_rag.database.models import Dataset, Experiment, Feature, Program, Signature
from amprenta_rag.ingestion.features.postgres_linking import link_features_to_postgres_items
from amprenta_rag.ingestion.postgres_program_experiment_linking import (
    link_dataset_to_programs_and_experiments_in_postgres,
)
from amprenta_rag.ingestion.postgres_signature_linking import link_signature_to_dataset_in_postgres
from scripts.dashboard.db_session import db_session


def render_management_page() -> None:
    """
    Render the Data Management page for editing metadata and linking entities.

    Features:
    - Edit dataset/program/experiment metadata
    - Link entities (datasets → programs, features → datasets)
    - Bulk operations
    """
    st.header("⚙️ Data Management")
    st.markdown("Edit metadata, link entities, and manage your data.")

    # Tab selection
    tab1, tab2, tab3 = st.tabs(["Link Entities", "Edit Metadata", "Bulk Operations"])

    # Tab 1: Link Entities
    with tab1:
        st.subheader("Link Entities")
        st.markdown("Create relationships between datasets, programs, experiments, and features.")

        link_type = st.selectbox(
            "Select Link Type",
            [
                "Dataset → Program",
                "Dataset → Experiment",
                "Feature → Dataset",
                "Signature → Dataset",
            ],
        )

        with db_session() as db:
            if link_type == "Dataset → Program":
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()
                programs = db.query(Program).order_by(Program.name).limit(200).all()

                if datasets and programs:
                    selected_dataset = st.selectbox(
                        "Select Dataset", [(d.id, d.name) for d in datasets], format_func=lambda x: x[1]
                    )
                    selected_programs = st.multiselect(
                        "Select Programs to Link", [(p.id, p.name) for p in programs], format_func=lambda x: x[1]
                    )

                    if st.button("🔗 Link Dataset to Programs", type="primary"):
                        if not selected_programs:
                            st.warning("Please select at least one program to link.")
                        else:
                            try:
                                with st.spinner("Linking dataset to programs..."):
                                    with db_session() as link_db:
                                        program_uuids = [p[0] for p in selected_programs]
                                        result = link_dataset_to_programs_and_experiments_in_postgres(
                                            dataset_id=selected_dataset[0],
                                            program_ids=program_uuids,
                                            db=link_db,
                                        )

                                        if result["programs_linked"] > 0:
                                            st.success(
                                                f"✅ Successfully linked dataset to {result['programs_linked']} program(s)!"
                                            )
                                            st.rerun()
                                        else:
                                            st.info("Dataset was already linked to all selected programs.")
                            except Exception as e:
                                st.error(f"❌ Error linking dataset to programs: {str(e)}")
                                st.exception(e)
                else:
                    if not datasets:
                        st.warning("No datasets available.")
                    if not programs:
                        st.warning("No programs available.")

            elif link_type == "Dataset → Experiment":
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()
                experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()

                if datasets and experiments:
                    selected_dataset = st.selectbox(
                        "Select Dataset",
                        [(d.id, d.name) for d in datasets],
                        format_func=lambda x: x[1],
                        key="ds_exp_ds",
                    )
                    selected_experiments = st.multiselect(
                        "Select Experiments to Link",
                        [(e.id, e.name) for e in experiments],
                        format_func=lambda x: x[1],
                        key="ds_exp_exp",
                    )

                    if st.button("🔗 Link Dataset to Experiments", type="primary", key="ds_exp_btn"):
                        if not selected_experiments:
                            st.warning("Please select at least one experiment to link.")
                        else:
                            try:
                                with st.spinner("Linking dataset to experiments..."):
                                    with db_session() as link_db:
                                        experiment_uuids = [e[0] for e in selected_experiments]
                                        result = link_dataset_to_programs_and_experiments_in_postgres(
                                            dataset_id=selected_dataset[0],
                                            experiment_ids=experiment_uuids,
                                            db=link_db,
                                        )

                                        if result["experiments_linked"] > 0:
                                            st.success(
                                                f"✅ Successfully linked dataset to {result['experiments_linked']} experiment(s)!"
                                            )
                                            st.rerun()
                                        else:
                                            st.info("Dataset was already linked to all selected experiments.")
                            except Exception as e:
                                st.error(f"❌ Error linking dataset to experiments: {str(e)}")
                                st.exception(e)
                else:
                    if not datasets:
                        st.warning("No datasets available.")
                    if not experiments:
                        st.warning("No experiments available.")

            elif link_type == "Feature → Dataset":
                features = db.query(Feature).order_by(Feature.name).limit(100).all()  # Limit for performance
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()

                if features and datasets:
                    feature_search = st.text_input("Search Feature", key="feat_search")
                    if feature_search:
                        features = [f for f in features if feature_search.lower() in f.name.lower()]

                    selected_feature = st.selectbox(
                        "Select Feature",
                        [(f.id, f.name) for f in features[:50]],  # Limit display
                        format_func=lambda x: x[1],
                        key="feat_ds_feat",
                    )
                    selected_datasets = st.multiselect(
                        "Select Datasets to Link",
                        [(d.id, d.name) for d in datasets],
                        format_func=lambda x: x[1],
                        key="feat_ds_ds",
                    )

                    if st.button("🔗 Link Feature to Datasets", type="primary", key="feat_ds_btn"):
                        if not selected_datasets:
                            st.warning("Please select at least one dataset to link.")
                        else:
                            try:
                                with st.spinner("Linking feature to datasets..."):
                                    feature = db.query(Feature).filter(Feature.id == selected_feature[0]).first()
                                    if feature:
                                        feature_name = feature.name

                                        linked_count = 0
                                        with db_session() as link_db:
                                            for dataset_tuple in selected_datasets:
                                                try:
                                                    link_features_to_postgres_items(
                                                        feature_names=[feature_name],
                                                        item_id=dataset_tuple[0],
                                                        item_type="dataset",
                                                        db=link_db,
                                                    )
                                                    linked_count += 1
                                                except Exception as e:
                                                    st.warning(f"Error linking to dataset {dataset_tuple[1]}: {str(e)}")

                                        if linked_count > 0:
                                            st.success(f"✅ Successfully linked feature to {linked_count} dataset(s)!")
                                            st.rerun()
                                        else:
                                            st.warning("No datasets were linked. Check for errors above.")
                                    else:
                                        st.error("Feature not found.")
                            except Exception as e:
                                st.error(f"❌ Error linking feature to datasets: {str(e)}")
                                st.exception(e)
                else:
                    if not features:
                        st.warning("No features available.")
                    if not datasets:
                        st.warning("No datasets available.")

            elif link_type == "Signature → Dataset":
                signatures = db.query(Signature).order_by(Signature.name).limit(200).all()
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()

                if signatures and datasets:
                    selected_signature = st.selectbox(
                        "Select Signature",
                        [(s.id, s.name) for s in signatures],
                        format_func=lambda x: x[1],
                        key="sig_ds_sig",
                    )
                    selected_datasets = st.multiselect(
                        "Select Datasets to Link",
                        [(d.id, d.name) for d in datasets],
                        format_func=lambda x: x[1],
                        key="sig_ds_ds",
                    )

                    if st.button("🔗 Link Signature to Datasets", type="primary", key="sig_ds_btn"):
                        if not selected_datasets:
                            st.warning("Please select at least one dataset to link.")
                        else:
                            try:
                                with st.spinner("Linking signature to datasets..."):
                                    linked_count = 0
                                    with db_session() as link_db:
                                        for dataset_tuple in selected_datasets:
                                            try:
                                                success = link_signature_to_dataset_in_postgres(
                                                    signature_id=selected_signature[0],
                                                    dataset_id=dataset_tuple[0],
                                                    db=link_db,
                                                )
                                                if success:
                                                    linked_count += 1
                                            except Exception as e:
                                                st.warning(f"Error linking to dataset {dataset_tuple[1]}: {str(e)}")

                                    if linked_count > 0:
                                        st.success(f"✅ Successfully linked signature to {linked_count} dataset(s)!")
                                        st.rerun()
                                    else:
                                        st.warning("No datasets were linked. Check for errors above.")
                            except Exception as e:
                                st.error(f"❌ Error linking signature to datasets: {str(e)}")
                                st.exception(e)
                else:
                    if not signatures:
                        st.warning("No signatures available.")
                    if not datasets:
                        st.warning("No datasets available.")

    # Tab 2: Edit Metadata
    with tab2:
        st.subheader("Edit Metadata")
        st.markdown("Update metadata for datasets, programs, experiments, and signatures.")

        entity_type = st.selectbox(
            "Select Entity Type", ["Dataset", "Program", "Experiment", "Signature"], key="edit_type"
        )

        with db_session() as db:
            if entity_type == "Dataset":
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()
                if datasets:
                    selected = st.selectbox(
                        "Select Dataset", [(d.id, d.name) for d in datasets], format_func=lambda x: x[1], key="edit_ds"
                    )
                    dataset_id = selected[0]
                    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()

                    if dataset:
                        st.markdown("### Current Metadata")
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**Name:** {dataset.name}")
                            st.write(f"**Omics Type:** {dataset.omics_type}")
                            st.write(f"**ID:** `{dataset.id}`")
                        with col2:
                            st.write(f"**Created:** {dataset.created_at}")
                            st.write(f"**Updated:** {dataset.updated_at}")

                        st.markdown("### Edit Metadata")
                        new_name = st.text_input("Name", value=dataset.name, key="edit_ds_name")
                        new_description = st.text_area(
                            "Description", value=dataset.description or "", key="edit_ds_desc"
                        )
                        new_disease = st.text_input(
                            "Disease (comma-separated)", value=", ".join(dataset.disease or []), key="edit_ds_disease"
                        )
                        new_sample_type = st.text_input(
                            "Sample Type (comma-separated)", value=", ".join(dataset.sample_type or []), key="edit_ds_sample"
                        )
                        new_organism = st.text_input(
                            "Organism (comma-separated)", value=", ".join(dataset.organism or []), key="edit_ds_org"
                        )
                        new_omics = st.selectbox(
                            "Omics Type",
                            ["lipidomics", "metabolomics", "proteomics", "transcriptomics", "other"],
                            index=["lipidomics", "metabolomics", "proteomics", "transcriptomics", "other"].index(
                                (dataset.omics_type or "other")
                            )
                            if (dataset.omics_type or "other") in ["lipidomics", "metabolomics", "proteomics", "transcriptomics", "other"]
                            else 4,
                            key="edit_ds_omics",
                        )
                        new_methods = st.text_area("Methods", value=dataset.methods or "", key="edit_ds_methods")
                        new_summary = st.text_area("Summary", value=dataset.summary or "", key="edit_ds_summary")
                        new_results = st.text_area("Results", value=dataset.results or "", key="edit_ds_results")
                        new_conclusions = st.text_area(
                            "Conclusions", value=dataset.conclusions or "", key="edit_ds_conclusions"
                        )
                        new_data_origin = st.text_input(
                            "Data Origin", value=dataset.data_origin or "", key="edit_ds_origin"
                        )
                        new_source_type = st.text_input(
                            "Dataset Source Type", value=dataset.dataset_source_type or "", key="edit_ds_source_type"
                        )

                        def _to_list(val: str):
                            return [x.strip() for x in val.split(",") if x.strip()] if val else []

                        if st.button("💾 Save Changes", type="primary", key="edit_ds_save"):
                            update_data = DatasetUpdate(
                                name=new_name if new_name != dataset.name else None,
                                description=(
                                    new_description if new_description != (dataset.description or "") else None
                                ),
                                disease=_to_list(new_disease) if _to_list(new_disease) != (dataset.disease or []) else None,
                                sample_type=_to_list(new_sample_type)
                                if _to_list(new_sample_type) != (dataset.sample_type or [])
                                else None,
                                organism=_to_list(new_organism)
                                if _to_list(new_organism) != (dataset.organism or [])
                                else None,
                                omics_type=new_omics if new_omics != (dataset.omics_type or "") else None,
                                methods=new_methods if new_methods != (dataset.methods or "") else None,
                                summary=new_summary if new_summary != (dataset.summary or "") else None,
                                results=new_results if new_results != (dataset.results or "") else None,
                                conclusions=new_conclusions
                                if new_conclusions != (dataset.conclusions or "")
                                else None,
                                data_origin=new_data_origin
                                if new_data_origin != (dataset.data_origin or "")
                                else None,
                                dataset_source_type=new_source_type
                                if new_source_type != (dataset.dataset_source_type or "")
                                else None,
                            )

                            if not any(
                                [
                                    update_data.name,
                                    update_data.description,
                                    update_data.disease,
                                    update_data.sample_type,
                                    update_data.organism,
                                    update_data.omics_type,
                                    update_data.methods,
                                    update_data.summary,
                                    update_data.results,
                                    update_data.conclusions,
                                    update_data.data_origin,
                                    update_data.dataset_source_type,
                                ]
                            ):
                                st.info("No changes detected.")
                            else:
                                try:
                                    with st.spinner("Saving changes..."):
                                        with db_session() as update_db:
                                            updated = update_dataset(
                                                db=update_db,
                                                dataset_id=dataset_id,
                                                dataset=update_data,
                                            )

                                            if updated:
                                                st.success("✅ Dataset metadata updated successfully!")
                                                st.rerun()
                                            else:
                                                st.error("Failed to update dataset.")
                                except Exception as e:
                                    st.error(f"❌ Error updating dataset: {str(e)}")
                                    st.exception(e)

                        st.markdown("### Raw Metadata")
                        with st.expander("View raw metadata", expanded=False):
                            st.json(
                                {
                                    "id": str(dataset.id),
                                    "name": dataset.name,
                                    "description": dataset.description,
                                    "omics_type": dataset.omics_type,
                                    "disease": dataset.disease,
                                    "sample_type": dataset.sample_type,
                                    "organism": dataset.organism,
                                    "file_paths": dataset.file_paths,
                                    "external_ids": dataset.external_ids,
                                    "methods": dataset.methods,
                                    "summary": dataset.summary,
                                    "results": dataset.results,
                                    "conclusions": dataset.conclusions,
                                    "data_origin": dataset.data_origin,
                                    "dataset_source_type": dataset.dataset_source_type,
                                    "created_at": str(dataset.created_at),
                                    "updated_at": str(dataset.updated_at),
                                }
                            )
                else:
                    st.warning("No datasets available.")

            elif entity_type == "Program":
                programs = db.query(Program).order_by(Program.name).limit(200).all()
                if programs:
                    selected = st.selectbox(
                        "Select Program",
                        [(p.id, p.name) for p in programs],
                        format_func=lambda x: x[1],
                        key="edit_prog",
                    )
                    program_id = selected[0]
                    program = db.query(Program).filter(Program.id == program_id).first()

                    if program:
                        st.markdown("### Current Metadata")
                        st.write(f"**Name:** {program.name}")
                        st.write(f"**Description:** {program.description or 'N/A'}")
                        st.write(f"**ID:** `{program.id}`")

                        st.markdown("### Edit Metadata")
                        new_name = st.text_input("Name", value=program.name, key="edit_prog_name")
                        new_description = st.text_area(
                            "Description", value=program.description or "", key="edit_prog_desc"
                        )

                        if st.button("💾 Save Changes", type="primary", key="edit_prog_save"):
                            if new_name == program.name and new_description == (program.description or ""):
                                st.info("No changes detected.")
                            else:
                                try:
                                    with st.spinner("Saving changes..."):
                                        with db_session() as update_db:
                                            update_data = ProgramUpdate(
                                                name=new_name if new_name != program.name else None,
                                                description=(
                                                    new_description
                                                    if new_description != (program.description or "")
                                                    else None
                                                ),
                                            )
                                            updated = update_program(
                                                db=update_db,
                                                program_id=program_id,
                                                program=update_data,
                                            )

                                            if updated:
                                                st.success("✅ Program metadata updated successfully!")
                                                st.rerun()
                                            else:
                                                st.error("Failed to update program.")
                                except Exception as e:
                                    st.error(f"❌ Error updating program: {str(e)}")
                                    st.exception(e)
                else:
                    st.warning("No programs available.")

            elif entity_type == "Experiment":
                experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()
                if experiments:
                    selected = st.selectbox(
                        "Select Experiment",
                        [(e.id, e.name) for e in experiments],
                        format_func=lambda x: x[1],
                        key="edit_exp",
                    )
                    experiment_id = selected[0]
                    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

                    if experiment:
                        st.markdown("### Current Metadata")
                        st.write(f"**Name:** {experiment.name}")
                        st.write(f"**Type:** {experiment.type or 'N/A'}")
                        st.write(f"**Description:** {experiment.description or 'N/A'}")
                        st.write(f"**ID:** `{experiment.id}`")

                        st.markdown("### Edit Metadata")
                        new_name = st.text_input("Name", value=experiment.name, key="edit_exp_name")
                        new_type = st.text_input("Type", value=experiment.type or "", key="edit_exp_type")
                        new_description = st.text_area(
                            "Description", value=experiment.description or "", key="edit_exp_desc"
                        )

                        if st.button("💾 Save Changes", type="primary", key="edit_exp_save"):
                            if (
                                new_name == experiment.name
                                and new_type == (experiment.type or "")
                                and new_description == (experiment.description or "")
                            ):
                                st.info("No changes detected.")
                            else:
                                try:
                                    with st.spinner("Saving changes..."):
                                        with db_session() as update_db:
                                            update_data = ExperimentUpdate(
                                                name=new_name if new_name != experiment.name else None,
                                                type=new_type if new_type != (experiment.type or "") else None,
                                                description=(
                                                    new_description
                                                    if new_description != (experiment.description or "")
                                                    else None
                                                ),
                                            )
                                            updated = update_experiment(
                                                db=update_db,
                                                experiment_id=experiment_id,
                                                experiment=update_data,
                                            )

                                            if updated:
                                                st.success("✅ Experiment metadata updated successfully!")
                                                st.rerun()
                                            else:
                                                st.error("Failed to update experiment.")
                                except Exception as e:
                                    st.error(f"❌ Error updating experiment: {str(e)}")
                                    st.exception(e)
                else:
                    st.warning("No experiments available.")

            else:  # Signature
                st.info("Signature editing is coming soon.")

    # Tab 3: Bulk Operations
    with tab3:
        st.subheader("Bulk Operations")
        st.markdown("Perform operations on multiple entities at once.")

        operation_type = st.selectbox(
            "Select Operation",
            [
                "Bulk Link Datasets to Program",
                "Bulk Export",
                "Bulk Delete (Coming Soon)",
            ],
        )

        if operation_type == "Bulk Link Datasets to Program":
            with db_session() as db:
                datasets = db.query(Dataset).order_by(Dataset.name).limit(200).all()
                programs = db.query(Program).order_by(Program.name).limit(200).all()

                if datasets and programs:
                    selected_datasets = st.multiselect(
                        "Select Datasets", [(d.id, d.name) for d in datasets], format_func=lambda x: x[1], key="bulk_ds"
                    )
                    selected_program = st.selectbox(
                        "Select Program",
                        [(p.id, p.name) for p in programs],
                        format_func=lambda x: x[1],
                        key="bulk_prog",
                    )

                    if st.button("🔗 Bulk Link", type="primary", key="bulk_link"):
                        if not selected_datasets:
                            st.warning("Please select at least one dataset to link.")
                        elif not selected_program:
                            st.warning("Please select a program to link to.")
                        else:
                            try:
                                with st.spinner(f"Linking {len(selected_datasets)} dataset(s) to program..."):
                                    program_id = selected_program[0]
                                    linked_count = 0
                                    errors = []

                                    with db_session() as link_db:
                                        for dataset_tuple in selected_datasets:
                                            try:
                                                result = link_dataset_to_programs_and_experiments_in_postgres(
                                                    dataset_id=dataset_tuple[0],
                                                    program_ids=[program_id],
                                                    db=link_db,
                                                )
                                                if result["programs_linked"] > 0:
                                                    linked_count += 1
                                            except Exception as e:
                                                errors.append(f"{dataset_tuple[1]}: {str(e)}")

                                    if linked_count > 0:
                                        st.success(
                                            f"✅ Successfully linked {linked_count}/{len(selected_datasets)} dataset(s) to program!"
                                        )
                                        if errors:
                                            st.warning(f"Errors for {len(errors)} dataset(s):")
                                            for error in errors[:5]:  # Show first 5 errors
                                                st.text(error)
                                        st.rerun()
                                    else:
                                        st.warning(
                                            "No datasets were linked. They may already be linked to this program."
                                        )
                            except Exception as e:
                                st.error(f"❌ Error during bulk linking: {str(e)}")
                                st.exception(e)
                else:
                    if not datasets:
                        st.warning("No datasets available.")
                    if not programs:
                        st.warning("No programs available.")

        elif operation_type == "Bulk Export":
            st.info("Bulk export functionality is available on individual entity pages (CSV download buttons).")

        else:
            st.warning("This operation is not yet available in the UI.")
