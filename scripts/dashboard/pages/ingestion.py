"""Data Ingestion page for the Streamlit dashboard."""

from __future__ import annotations

import tempfile
from pathlib import Path

import streamlit as st

from amprenta_rag.database.models import Experiment, Program
from amprenta_rag.domain.omics import OmicsDatasetIngestRequest
from amprenta_rag.ingestion.local_documents import ingest_local_document
from amprenta_rag.ingestion.omics_service import ingest_dataset_from_file
from scripts.dashboard.auth import require_admin, require_auth
from scripts.dashboard.db_session import db_session
from scripts.dashboard.utils.errors import render_ingest_error
from scripts.dashboard.utils.ingest_summary import render_dataset_ingest_summary


def error_message(message: str, what_next: str = None):
    st.error(message)
    if what_next:
        st.info(f"What to do next: {what_next}")


def render_ingestion_page() -> None:
    """
    Render the Data Ingestion page for uploading and ingesting datasets.

    Features:
    - File upload (CSV, TSV, mwTab)
    - Omics type selection
    - Metadata input form
    - Program/Experiment linking
    - Progress tracking
    """
    user = require_auth()
    require_admin(user)
    st.header("📥 Data Ingestion")
    st.markdown("Upload and ingest datasets into the platform.")

    # Tab selection
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(
        [
            "Dataset Upload",
            "Signature Upload",
            "Repository Harvest",
            "Chemistry Data",
            "Email/Zotero",
            "Local Documents",
        ]
    )

    # Tab 1: Dataset Upload
    with tab1:
        st.subheader("Upload Dataset File")
        uploaded_file = st.file_uploader(
            "Upload dataset file (CSV/TSV)", type=["csv", "tsv", "txt"], help="Upload dataset file (required)"
        )

        # Load program and experiment names (extract while session is active)
        with db_session() as db:
            programs = db.query(Program).order_by(Program.name).limit(200).all()
            experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()
        with db_session() as db:
            programs = db.query(Program).order_by(Program.name).limit(200).all()
            experiments = db.query(Experiment).order_by(Experiment.name).limit(200).all()
            program_names = [p.name for p in programs]
            experiment_names = [e.name for e in experiments]

        with st.form("dataset_upload_form", clear_on_submit=False):
            dataset_name = st.text_input("Dataset name*", help="Required", key="dataset_name")
            omics_type = st.selectbox(
                "Omics type*",
                ["lipidomics", "metabolomics", "proteomics", "transcriptomics", "multi-omics"],
                key="omics_type",
            )
            program_opt = st.selectbox("Program", ["(None)"] + program_names, key="program")
            experiment_opt = st.selectbox("Experiment", ["(None)"] + experiment_names, key="experiment")

            # Disease/Condition
            disease_choices = ["ALS", "AD", "Control", "PD", "Cancer", "Other"]
            diseases = st.multiselect("Disease / Condition*", disease_choices, help="Can select multiple")
            species = st.selectbox("Species*", ["Human", "Mouse", "Rat", "Other"], key="species")
            matrix = st.selectbox(
                "Matrix / Sample type*", ["CSF", "Plasma", "Serum", "Tissue", "Cell line", "Other"], key="matrix"
            )
            cohort_label = st.text_input("Cohort / Group label (optional)")
            design_summary = st.text_area("Design summary* (1–3 sentences)", height=100)
            normalization_status = st.selectbox(
                "Normalization status*",
                ["raw", "normalized", "log-transformed", "batch-corrected", "unknown"],
                key="norm_status",
            )
            value_type = st.selectbox(
                "Value type*",
                ["intensity", "counts", "log2 fold-change", "z-score", "p-value/FDR", "other"],
                key="value_type",
            )
            reference_group = st.text_input("Reference/baseline group (optional)", key="ref_group")

            # Omics-specific
            tech_details = {}
            if omics_type == "transcriptomics":
                tech_details["platform"] = st.text_input("Transcriptomics platform")
                tech_details["genome_build"] = st.text_input("Annotation/genome build")
            elif omics_type == "proteomics":
                tech_details["quantitation"] = st.selectbox(
                    "Quantitation", ["(None)", "LFQ", "TMT", "DIA", "Other"], key="quant"
                )
                tech_details["pipeline"] = st.text_input("Pipeline (e.g. MaxQuant)")
            elif omics_type in ("lipidomics", "metabolomics"):
                tech_details["acquisition_mode"] = st.text_input("Acquisition mode")
                tech_details["ionization_mode"] = st.text_input("Ionization mode")
                tech_details["internal_standard"] = st.text_input("Internal standard")

            include_discovery = st.checkbox("Include in automated discovery/validation", value=True)
            confidential = st.checkbox("Confidential / internal-only", value=False)

            submit = st.form_submit_button("🚀 Ingest Dataset", type="primary")

        if submit:
            if not all(
                [
                    dataset_name,
                    omics_type,
                    diseases,
                    species,
                    matrix,
                    design_summary,
                    normalization_status,
                    value_type,
                    uploaded_file,
                ]
            ):
                st.error("Please complete all required fields and upload a file.")
                st.stop()
            uploads_dir = Path("data/uploads")
            uploads_dir.mkdir(parents=True, exist_ok=True)
            file_path = uploads_dir / uploaded_file.name
            progress = st.progress(0, text="Saving uploaded file...")
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getvalue())
            progress.progress(0.2, text="Preparing ingestion request...")

            # Resolve IDs
            program_id, experiment_id = None, None
            for p in programs:
                if p.name == program_opt:
                    program_id = str(p.id)
            for e in experiments:
                if e.name == experiment_opt:
                    experiment_id = str(e.id)

            # Build semantic metadata
            semantic_metadata = dict(
                diseases=diseases,
                species=species,
                matrix=matrix,
                cohort_label=cohort_label,
                design_summary=design_summary,
                normalization_status=normalization_status,
                value_type=value_type,
                reference_group=reference_group,
                omics_type=omics_type,
                include_in_discovery=include_discovery,
                confidential=confidential,
                technical_details=tech_details,
            )

            req = OmicsDatasetIngestRequest(
                omics_type=omics_type,
                dataset_path=str(file_path),
                metadata=semantic_metadata,
                extra_args={
                    "name": dataset_name,
                    "description": design_summary,
                    "program_id": program_id,
                    "experiment_id": experiment_id,
                },
            )
            try:
                progress.progress(0.4, text="Ingesting dataset...")
                result = ingest_dataset_from_file(req)
                progress.progress(0.7, text="Refreshing status...")
                # Fetch new dataset for summary details (ID, QC, etc.)
                from amprenta_rag.database.models import Dataset

                with db_session() as db:
                    ds = db.query(Dataset).filter_by(id=result).first()
                progress.progress(1.0, text="Done!")
                render_dataset_ingest_summary(ds)
                st.balloons()
            except Exception as e:
                render_ingest_error(e)
                progress.empty()

    # Tab 2: Signature Upload
    with tab2:
        st.subheader("Upload Signature File")
        st.markdown("Upload a signature definition file (TSV/CSV format).")

        uploaded_sig = st.file_uploader(
            "Choose a signature file",
            type=["tsv", "csv", "txt"],
            key="signature_upload",
            help="Upload TSV or CSV format signature files",
        )

        if uploaded_sig:
            st.info(f"📄 File: {uploaded_sig.name} ({uploaded_sig.size:,} bytes)")

            # Preview file
            if st.checkbox("Preview file contents"):
                try:
                    content = uploaded_sig.read().decode("utf-8")
                    st.text_area("File Preview", content, height=200)
                    uploaded_sig.seek(0)  # Reset file pointer
                except Exception as e:
                    st.error(f"Error reading file: {e}")

            # Signature metadata
            col1, col2 = st.columns(2)
            with col1:
                signature_type = st.selectbox(
                    "Signature Type", ["Consortium", "Literature-derived", "Open Dataset", "Other"], key="sig_type"
                )
                data_ownership = st.text_input("Data Ownership", value="Public", key="sig_ownership")
            with col2:
                version = st.text_input("Version (optional)", key="sig_version")
                description = st.text_area("Description (optional)", height=100, key="sig_description")

            # Optional metadata
            disease_context = st.text_input(
                "Disease Context (optional, comma-separated)", help="Enter disease contexts separated by commas"
            )
            matrix = st.text_input("Matrix (optional, comma-separated)", help="Enter matrices separated by commas")

            if st.button("🚀 Ingest Signature", type="primary", key="ingest_sig"):
                if not uploaded_sig:
                    st.error("Please upload a signature file first.")
                else:
                    with st.spinner("Ingesting signature... This may take a few minutes."):
                        try:
                            from amprenta_rag.ingestion.signature_ingestion import ingest_signature_from_file

                            # Save uploaded file to temporary location
                            with tempfile.NamedTemporaryFile(
                                delete=False, suffix=Path(uploaded_sig.name).suffix, mode="wb"
                            ) as tmp_file:
                                tmp_file.write(uploaded_sig.getvalue())
                                tmp_path = tmp_file.name

                            # Parse optional fields
                            disease_list = [d.strip() for d in disease_context.split(",")] if disease_context else None
                            matrix_list = [m.strip() for m in matrix.split(",")] if matrix else None

                            result = ingest_signature_from_file(
                                signature_path=Path(tmp_path),
                                signature_type=signature_type,
                                data_ownership=data_ownership,
                                version=version if version else None,
                                description=description if description else None,
                                disease_context=disease_list,
                                matrix=matrix_list,
                            )

                            # Clean up temp file
                            Path(tmp_path).unlink()

                            st.success("✅ Signature ingested successfully!")
                            st.info(f"Signature Page ID: `{result.get('signature_page_id', 'N/A')}`")
                            if result.get("components_created"):
                                st.info(f"Components Created: {result['components_created']}")
                            st.balloons()

                        except Exception as e:
                            st.error(f"❌ Signature ingestion failed: {str(e)}")
                            st.exception(e)

    # Tab 3: Repository Harvest
    with tab3:
        st.subheader("Harvest Repository Study")
        st.markdown("Harvest studies from public repositories (GEO, PRIDE, MetaboLights, Metabolomics Workbench).")

        repository = st.selectbox(
            "Select Repository",
            ["GEO", "PRIDE", "MetaboLights", "Metabolomics Workbench"],
            help="Select the repository to harvest from",
        )

        study_id = st.text_input(
            "Study ID", help=f"Enter the {repository} study identifier (e.g., GSE12345, PXD000001, MTBLS123, ST000123)"
        )

        col1, col2 = st.columns(2)
        with col1:
            create_postgres = st.checkbox(
                "Create Postgres Dataset", value=True, help="Create dataset in Postgres database"
            )
        with col2:
            ingest_after = st.checkbox(
                "Ingest After Creation", value=False, help="Automatically ingest the dataset after creating it"
            )

        if st.button("🌾 Harvest Study", type="primary"):
            if not study_id:
                st.error("Please enter a study ID.")
            else:
                with st.spinner(f"Harvesting {repository} study {study_id}... This may take a few minutes."):
                    try:
                        from scripts.harvest_repository_study import harvest_study

                        dataset_id, notion_page_id = harvest_study(
                            study_id=study_id,
                            repository=repository,
                            create_notion=False,  # Postgres-first
                            create_postgres=create_postgres,
                            ingest=ingest_after,
                            dry_run=False,
                        )

                        if dataset_id:
                            st.success("✅ Study harvested successfully!")
                            st.info(f"Dataset ID: `{dataset_id}`")
                            if notion_page_id:
                                st.info(f"Notion Page ID: `{notion_page_id}`")
                            st.balloons()
                        else:
                            st.warning("⚠️ Harvest completed but no dataset ID returned. Check logs for details.")

                    except Exception as e:
                        st.error(f"❌ Harvest failed: {str(e)}")
                        st.exception(e)

        st.info(
            'To search repositories by disease, species, or advanced metadata, use the “Repositories” page in the sidebar navigation.'
        )

    # Tab 5: Email/Zotero Ingestion
    with tab5:
        st.subheader("Email & Zotero Ingestion")

        ingestion_type = st.radio("Select Ingestion Type", ["Gmail", "Zotero Collection"], horizontal=True)

        if ingestion_type == "Gmail":
            st.markdown("### Gmail Ingestion")
            st.info("Gmail ingestion requires OAuth setup. Use the CLI script for initial setup.")

            col1, col2 = st.columns(2)
            with col1:
                max_emails = st.number_input(
                    "Maximum Emails", min_value=1, max_value=1000, value=50, help="Maximum number of emails to ingest"
                )
            with col2:
                query = st.text_input(
                    "Gmail Query (optional)", help="Gmail search query (e.g., 'from:example@gmail.com')"
                )

            if st.button("📧 Ingest Gmail", type="primary"):
                st.info("Gmail ingestion via UI is coming soon. Please use `scripts/ingest_gmail.py` for now.")
                # TODO: Implement Gmail ingestion UI

        else:  # Zotero
            st.markdown("### Zotero Collection Ingestion")

            with db_session() as db:
                # Get available collections (would need Zotero API integration)
                collection_id = st.text_input("Zotero Collection ID", help="Enter the Zotero collection ID to ingest")

                if collection_id:
                    if st.button("📚 Ingest Zotero Collection", type="primary"):
                        with st.spinner("Ingesting Zotero collection... This may take a few minutes."):
                            try:

                                # Note: This would need to be adapted for UI use
                                st.info(
                                    "Zotero ingestion via UI is coming soon. Please use `scripts/ingest_zotero_postgres.py` for now."
                                )
                                # TODO: Implement Zotero ingestion UI

                            except Exception as e:
                                st.error(f"❌ Zotero ingestion failed: {str(e)}")
                                st.exception(e)

    # Tab 4: Chemistry Data Upload
    with tab4:
        st.subheader("Upload Chemistry Data")
        st.markdown("Upload compounds, HTS campaigns, and biochemical assay results.")

        chemistry_type = st.radio(
            "Select Data Type", ["HTS Campaign", "Biochemical Results", "Compound List"], horizontal=True
        )

        if chemistry_type == "HTS Campaign":
            st.markdown("### HTS Campaign Upload")

            # Campaign metadata file
            metadata_file = st.file_uploader(
                "Upload Campaign Metadata (CSV/TSV)",
                type=["csv", "tsv"],
                key="hts_metadata",
                help="CSV/TSV file with campaign metadata (campaign_name, description, assay_type, target, etc.)",
            )

            if metadata_file:
                st.info(f"📄 Metadata file: {metadata_file.name}")

                # Hit list file
                hit_list_file = st.file_uploader(
                    "Upload Hit List (CSV/TSV)",
                    type=["csv", "tsv"],
                    key="hts_hitlist",
                    help="CSV/TSV file with SMILES and assay values",
                )

                if hit_list_file:
                    st.info(f"📄 Hit list file: {hit_list_file.name}")

                    # Configuration
                    col1, col2 = st.columns(2)
                    with col1:
                        smiles_column = st.text_input(
                            "SMILES Column Name", value="SMILES", help="Name of the column containing SMILES strings"
                        )
                        value_column = st.text_input(
                            "Value Column Name", value="value", help="Name of the column containing assay values"
                        )
                    with col2:
                        vendor_id_column = st.text_input(
                            "Vendor ID Column (optional)", value="", help="Optional column name for vendor compound IDs"
                        )
                        hit_threshold = st.number_input(
                            "Hit Threshold (optional)",
                            min_value=0.0,
                            value=0.0,
                            help="Threshold for hit classification",
                        )

                    if st.button("🚀 Ingest HTS Campaign", type="primary"):
                        with st.spinner("Ingesting HTS campaign... This may take a few minutes."):
                            try:
                                from pathlib import Path

                                from amprenta_rag.ingestion.screening_ingestion import (
                                    ingest_hts_campaign_metadata,
                                    ingest_hts_hit_list,
                                )

                                # Save files temporarily
                                with tempfile.NamedTemporaryFile(
                                    delete=False, suffix=Path(metadata_file.name).suffix, mode="wb"
                                ) as tmp_meta:
                                    tmp_meta.write(metadata_file.getvalue())
                                    tmp_meta_path = tmp_meta.name

                                with tempfile.NamedTemporaryFile(
                                    delete=False, suffix=Path(hit_list_file.name).suffix, mode="wb"
                                ) as tmp_hits:
                                    tmp_hits.write(hit_list_file.getvalue())
                                    tmp_hits_path = tmp_hits.name

                                # Ingest metadata
                                campaign = ingest_hts_campaign_metadata(
                                    metadata_file=Path(tmp_meta_path),
                                )

                                # Ingest hit list
                                compounds_count = ingest_hts_hit_list(
                                    hit_list_file=Path(tmp_hits_path),
                                    campaign_id=campaign.campaign_id,
                                    smiles_column=smiles_column,
                                    vendor_id_column=vendor_id_column if vendor_id_column else None,
                                    value_column=value_column,
                                    hit_threshold=hit_threshold if hit_threshold > 0 else None,
                                )

                                # Clean up
                                Path(tmp_meta_path).unlink()
                                Path(tmp_hits_path).unlink()

                                st.success("✅ HTS campaign ingested successfully!")
                                st.info(f"Campaign ID: `{campaign.campaign_id}`")
                                st.info(f"Compounds ingested: {compounds_count}")
                                st.balloons()

                            except Exception as e:
                                st.error(f"❌ HTS campaign ingestion failed: {str(e)}")
                                st.exception(e)

        elif chemistry_type == "Biochemical Results":
            st.markdown("### Biochemical Assay Results Upload")

            results_file = st.file_uploader(
                "Upload Biochemical Results (CSV/TSV)",
                type=["csv", "tsv"],
                key="bio_results",
                help="CSV/TSV file with assay results (result_id, assay_name, target, IC50, EC50, etc.)",
            )

            if results_file:
                st.info(f"📄 Results file: {results_file.name}")

                # Preview
                if st.checkbox("Preview file contents"):
                    try:
                        import pandas as pd

                        df = pd.read_csv(results_file, sep="\t" if results_file.name.endswith(".tsv") else ",")
                        st.dataframe(df.head(10))
                        results_file.seek(0)
                    except Exception as e:
                        st.error(f"Error reading file: {e}")

                # Configuration
                col1, col2 = st.columns(2)
                with col1:
                    smiles_column = st.text_input(
                        "SMILES Column Name",
                        value="SMILES",
                        help="Name of the column containing SMILES strings",
                        key="bio_smiles_col",
                    )
                    assay_name = st.text_input(
                        "Assay Name",
                        value=Path(results_file.name).stem if results_file else "",
                        help="Name for this assay (defaults to filename)",
                        key="bio_assay_name",
                    )
                with col2:
                    target_column = st.text_input(
                        "Target Column (optional)",
                        value="target",
                        help="Name of the column containing target information",
                        key="bio_target_col",
                    )
                    units = st.text_input(
                        "Units", value="nM", help="Units for activity values (e.g., nM, μM)", key="bio_units"
                    )

                # Activity value columns
                st.markdown("#### Activity Value Columns")
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    ic50_column = st.text_input("IC50 Column", value="IC50", key="bio_ic50")
                with col2:
                    ec50_column = st.text_input("EC50 Column", value="EC50", key="bio_ec50")
                with col3:
                    ki_column = st.text_input("Ki Column", value="Ki", key="bio_ki")
                with col4:
                    kd_column = st.text_input("Kd Column", value="Kd", key="bio_kd")

                if st.button("🚀 Ingest Biochemical Results", type="primary"):
                    with st.spinner("Ingesting biochemical results... This may take a few minutes."):
                        try:
                            from pathlib import Path

                            from amprenta_rag.ingestion.screening_ingestion import ingest_biochemical_results

                            # Save file temporarily
                            with tempfile.NamedTemporaryFile(
                                delete=False, suffix=Path(results_file.name).suffix, mode="wb"
                            ) as tmp_file:
                                tmp_file.write(results_file.getvalue())
                                tmp_path = tmp_file.name

                            # Ingest results
                            results_count = ingest_biochemical_results(
                                results_file=Path(tmp_path),
                                smiles_column=smiles_column,
                                assay_name=assay_name if assay_name else None,
                            )

                            # Clean up
                            Path(tmp_path).unlink()

                            st.success("✅ Biochemical results ingested successfully!")
                            st.info(f"Results ingested: {results_count}")
                            st.balloons()

                        except Exception as e:
                            st.error(f"❌ Biochemical results ingestion failed: {str(e)}")
                            st.exception(e)

        else:  # Compound List
            st.markdown("### Compound List Upload")

            compound_file = st.file_uploader(
                "Upload Compound List (CSV/TSV)",
                type=["csv", "tsv"],
                key="compound_list",
                help="CSV/TSV file with compounds (SMILES, compound_id, etc.)",
            )

            if compound_file:
                st.info(f"📄 Compound file: {compound_file.name}")

                # Preview
                if st.checkbox("Preview file contents"):
                    try:
                        import pandas as pd

                        df = pd.read_csv(compound_file, sep="\t" if compound_file.name.endswith(".tsv") else ",")
                        st.dataframe(df.head(10))
                        compound_file.seek(0)
                    except Exception as e:
                        st.error(f"Error reading file: {e}")

                # Configuration
                col1, col2 = st.columns(2)
                with col1:
                    smiles_column = st.text_input(
                        "SMILES Column Name",
                        value="SMILES",
                        help="Name of the column containing SMILES strings",
                        key="compound_smiles_col",
                    )
                    compound_id_column = st.text_input(
                        "Compound ID Column (optional)",
                        value="",
                        help="Optional column name for compound IDs (auto-generated if not provided)",
                        key="compound_id_col",
                    )
                with col2:
                    vendor_id_column = st.text_input(
                        "Vendor ID Column (optional)",
                        value="",
                        help="Optional column name for vendor compound IDs",
                        key="compound_vendor_col",
                    )

                if st.button("🚀 Ingest Compounds", type="primary"):
                    with st.spinner("Ingesting compounds... This may take a few minutes."):
                        try:
                            import pandas as pd

                            from amprenta_rag.chemistry.registration import register_compound
                            from amprenta_rag.chemistry.normalization import (
                                compute_molecular_descriptors,
                                generate_compound_id,
                                normalize_smiles,
                            )

                            # Read compound file
                            df = pd.read_csv(compound_file, sep="\t" if compound_file.name.endswith(".tsv") else ",")

                            if smiles_column not in df.columns:
                                st.error(
                                    f"SMILES column '{smiles_column}' not found in file. Available columns: {', '.join(df.columns)}"
                                )
                            else:
                                compounds_created = 0
                                errors = []

                                progress_bar = st.progress(0)
                                total_rows = len(df)

                                for idx, row in df.iterrows():
                                    try:
                                        smiles = str(row[smiles_column]).strip()
                                        if not smiles or smiles == "nan":
                                            continue

                                        # Normalize compound
                                        canonical_smiles, inchi_key, molecular_formula = normalize_smiles(smiles)

                                        # Get compound ID
                                        if compound_id_column and compound_id_column in df.columns:
                                            compound_id = str(row[compound_id_column]).strip()
                                        else:
                                            compound_id = generate_compound_id(smiles)

                                        # Register compound (handles duplicates/upserts)
                                            registered_id = register_compound(
                                                name=compound_id,
                                                smiles=smiles,
                                                salt_form=None,
                                            )
                                        if registered_id:
                                            compounds_created += 1

                                    except Exception as e:
                                        errors.append(f"Row {idx + 1}: {str(e)}")

                                    # Update progress
                                    progress_bar.progress((idx + 1) / total_rows)

                                progress_bar.empty()

                                if errors:
                                    st.warning(f"⚠️ Completed with {len(errors)} error(s):")
                                    for error in errors[:10]:  # Show first 10 errors
                                        st.caption(f"  - {error}")
                                    if len(errors) > 10:
                                        st.caption(f"  ... and {len(errors) - 10} more errors")

                                st.success("✅ Compounds ingested successfully!")
                                st.info(f"Compounds processed: {compounds_created}")
                                if errors:
                                    st.info(f"Errors: {len(errors)}")
                                st.balloons()

                        except Exception as e:
                            st.error(f"❌ Compound ingestion failed: {str(e)}")
                            st.exception(e)

    with tab6:
        st.subheader("Local Document Ingestion")
        st.markdown("Upload one or more PDFs or text files for literature, protocols, or internal notes.")
        docs = st.file_uploader(
            "Upload documents",
            type=["pdf", "txt", "md"],
            accept_multiple_files=True,
            help="PDFs and plain text/markdown files.",
        )
        # Metadata form
        doc_type = st.selectbox("Document type", ["Literature", "Internal Note", "Protocol", "Other"])
        with db_session() as db:
            programs = db.query(Program).order_by(Program.name).all()
            # Extract program names inside session
            program_names = [p.name for p in programs]
        program_opt = st.selectbox("Program", ["(None)"] + program_names, key="doc_prog")
        disease_choices = ["ALS", "AD", "Control", "PD", "Cancer", "Other"]
        diseases = st.multiselect("Disease / Condition", disease_choices)
        tags = st.text_input("Tags (comma-separated)")
        year = st.text_input("Year (optional)")
        doi = st.text_input("DOI/PMID (optional)")
        journal = st.text_input("Journal (optional)")
        url = st.text_input("URL (optional)")
        submit_doc = st.button("Ingest Documents", type="primary")
        if submit_doc:
            if not docs:
                st.error("Please upload at least one document.")
                st.stop()
            upload_dir = Path("data/local_documents")
            upload_dir.mkdir(parents=True, exist_ok=True)
            for uploaded_file in docs:
                # Per-file title prompt
                title = st.text_input(f"Title for {uploaded_file.name} (optional)", value=Path(uploaded_file.name).stem)
                local_path = upload_dir / uploaded_file.name
                local_path.write_bytes(uploaded_file.getbuffer())
                # Resolve program_id
                program_id = None
                for p in programs:
                    if p.name == program_opt:
                        program_id = str(p.id)
                metadata = {
                    "program_id": program_id,
                    "diseases": diseases,
                    "tags": [t.strip() for t in tags.split(",") if t.strip()],
                    "year": year,
                    "doi": doi,
                    "journal": journal,
                    "url": url,
                }
                try:
                    doc_id = ingest_local_document(
                        file_path=local_path,
                        doc_type=doc_type.lower(),
                        title=title or Path(uploaded_file.name).stem,
                        metadata=metadata,
                    )
                    st.success(f"Document ingested: {title or uploaded_file.name} (id: {doc_id})")
                    st.info(f"💡 Go to **Literature** page to view document: `{doc_id}`")
                    st.code(f"Document ID: {doc_id}", language=None)
                except Exception as e:
                    from scripts.dashboard.utils.errors import render_ingest_error

                    render_ingest_error(e)
