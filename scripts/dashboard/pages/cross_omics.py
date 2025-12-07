"""Cross-Omics Summary page for the Streamlit dashboard."""

from __future__ import annotations

from uuid import UUID

import streamlit as st

from amprenta_rag.database.models import Dataset, Feature, Program, Signature
from amprenta_rag.query.cross_omics.dataset_summary_postgres import cross_omics_dataset_summary_postgres
from amprenta_rag.query.cross_omics.feature_summary_postgres import cross_omics_feature_summary_postgres
from amprenta_rag.query.cross_omics.program_summary_postgres import cross_omics_program_summary_postgres
from amprenta_rag.query.cross_omics.signature_summary_postgres import cross_omics_signature_summary_postgres
from scripts.dashboard.db_session import db_session


def render_cross_omics_page() -> None:
    """
    Render the Cross-Omics Summary page for generating multi-omics summaries.

    Features:
    - Program summary
    - Feature summary
    - Signature summary
    - Dataset summary
    """
    st.header("🔬 Cross-Omics Summary")
    st.markdown("Generate evidence-based summaries across multiple omics types.")

    # Summary type selection
    summary_type = st.radio(
        "Select summary type",
        ["Program", "Feature", "Signature", "Dataset"],
        horizontal=True,
    )

    if summary_type == "Program":
        st.subheader("Program Summary")

        with db_session() as db:
            # Get programs
            programs = db.query(Program).order_by(Program.name).all()

            if programs:
                program_options = {f"{p.name} ({p.id})": str(p.id) for p in programs}
                selected_program = st.selectbox(
                    "Select a program",
                    list(program_options.keys()),
                )
                program_id = UUID(program_options[selected_program])

                top_k = st.slider(
                    "Top K datasets per omics type", min_value=1, max_value=20, value=5, key="program_top_k"
                )

                if st.button("Generate Summary", type="primary"):
                    with st.spinner("Generating cross-omics summary..."):
                        try:
                            summary = cross_omics_program_summary_postgres(
                                program_id=program_id,
                                top_k_per_omics=top_k,
                            )
                            st.success("Summary generated!")
                            st.markdown("---")
                            st.markdown(summary)
                        except Exception as e:
                            st.error(f"Error generating summary: {str(e)}")
                            st.exception(e)
            else:
                st.info("No programs found. Create a program first.")

    elif summary_type == "Feature":
        st.subheader("Feature Summary")

        col1, col2 = st.columns(2)
        with col1:
            feature_name = st.text_input("Feature name (e.g., glucose, Cer(d18:1/16:0))")
        with col2:
            feature_type = st.selectbox(
                "Feature type",
                ["metabolite", "lipid", "gene", "protein"],
            )

        with db_session() as db:
            if feature_name:
                # Try to find feature
                feature = (
                    db.query(Feature)
                    .filter(
                        Feature.name.ilike(f"%{feature_name}%"),
                        Feature.feature_type == feature_type,
                    )
                    .first()
                )

                if feature:
                    st.info(f"Found feature: {feature.name} (ID: {feature.id})")
                    feature_id = feature.id
                else:
                    st.warning(f"Feature '{feature_name}' not found. Will search by name.")
                    feature_id = None
            else:
                feature_id = None

        col1, col2 = st.columns(2)
        with col1:
            top_k_datasets = st.slider("Top K datasets", min_value=1, max_value=50, value=20, key="feature_datasets")
        with col2:
            top_k_chunks = st.slider("Top K chunks", min_value=1, max_value=200, value=100, key="feature_chunks")

        if st.button("Generate Summary", type="primary"):
            if not feature_name:
                st.warning("Please enter a feature name.")
            else:
                with st.spinner("Generating cross-omics summary..."):
                    try:
                        summary = cross_omics_feature_summary_postgres(
                            feature_name=feature_name,
                            feature_type=feature_type,
                            feature_id=feature_id,
                            top_k_datasets=top_k_datasets,
                            top_k_chunks=top_k_chunks,
                        )
                        st.success("Summary generated!")
                        st.markdown("---")
                        st.markdown(summary)
                    except Exception as e:
                        st.error(f"Error generating summary: {str(e)}")
                        st.exception(e)

    elif summary_type == "Signature":
        st.subheader("Signature Summary")

        with db_session() as db:
            # Get signatures
            signatures = db.query(Signature).order_by(Signature.name).all()

            if signatures:
                signature_options = {f"{s.name} ({s.id})": str(s.id) for s in signatures}
                selected_signature = st.selectbox(
                    "Select a signature",
                    list(signature_options.keys()),
                )
                signature_id = UUID(signature_options[selected_signature])

                col1, col2 = st.columns(2)
                with col1:
                    top_k_datasets = st.slider(
                        "Top K datasets", min_value=1, max_value=50, value=20, key="sig_datasets"
                    )
                with col2:
                    top_k_chunks = st.slider("Top K chunks", min_value=1, max_value=200, value=100, key="sig_chunks")

                if st.button("Generate Summary", type="primary"):
                    with st.spinner("Generating cross-omics summary..."):
                        try:
                            summary = cross_omics_signature_summary_postgres(
                                signature_id=signature_id,
                                top_k_datasets=top_k_datasets,
                                top_k_chunks=top_k_chunks,
                            )
                            st.success("Summary generated!")
                            st.markdown("---")
                            st.markdown(summary)
                        except Exception as e:
                            st.error(f"Error generating summary: {str(e)}")
                            st.exception(e)
            else:
                st.info("No signatures found. Create a signature first.")

    elif summary_type == "Dataset":
        st.subheader("Dataset Summary")

        with db_session() as db:
            # Get datasets
            datasets = db.query(Dataset).order_by(Dataset.name).all()

            if datasets:
                dataset_options = {f"{d.name} ({d.id})": str(d.id) for d in datasets}
                selected_dataset = st.selectbox(
                    "Select a dataset",
                    list(dataset_options.keys()),
                )
                dataset_id = UUID(dataset_options[selected_dataset])

                top_k_chunks = st.slider("Top K chunks", min_value=1, max_value=200, value=100, key="dataset_chunks")

                if st.button("Generate Summary", type="primary"):
                    with st.spinner("Generating cross-omics summary..."):
                        try:
                            summary = cross_omics_dataset_summary_postgres(
                                dataset_id=dataset_id,
                                top_k_chunks=top_k_chunks,
                            )
                            st.success("Summary generated!")
                            st.markdown("---")
                            st.markdown(summary)
                        except Exception as e:
                            st.error(f"Error generating summary: {str(e)}")
                            st.exception(e)
            else:
                st.info("No datasets found. Ingest a dataset first.")
