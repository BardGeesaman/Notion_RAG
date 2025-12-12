"""Datasets page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from amprenta_rag.export.slide_generator import generate_dataset_slides
from scripts.dashboard.db_session import db_session


def get_qc_dict(ds):
    return ds.qc_metrics if hasattr(ds, "qc_metrics") and ds.qc_metrics else {}


def render_datasets_page() -> None:
    """
    Render the Datasets page with filtering, search, and detailed views.

    Features:
    - Filter by omics type
    - Search by name
    - Summary and detailed dataset views
    - Export functionality (CSV)
    """
    st.header("📁 Datasets")

    with db_session() as db:
        # Filters
        col1, col2, col3 = st.columns(3)
        with col1:
            # Get distinct omics types
            omics_types = [
                row[0]
                for row in db.query(Dataset.omics_type).distinct().limit(200).all()
                if row[0]
            ]
            omics_filter = st.selectbox(
                "Filter by Omics Type",
                ["All"] + omics_types,
            )
        with col2:
            # Use ingestion_status instead of qc_status (which doesn't exist)
            status_opts = [
                x[0]
                for x in db.query(Dataset.ingestion_status).distinct().limit(200).all()
                if x[0]
            ]
            qc_filter = st.selectbox("Ingestion Status", ["All"] + status_opts)
        with col3:
            search_term = st.text_input("Search by name", "")

        # Query datasets
        query = db.query(Dataset)

        if omics_filter != "All":
            query = query.filter(Dataset.omics_type == omics_filter)

        if qc_filter != "All":
            query = query.filter(Dataset.ingestion_status == qc_filter)

        if search_term:
            query = query.filter(Dataset.name.ilike(f"%{search_term}%"))

        datasets = query.order_by(Dataset.created_at.desc()).limit(200).all()

        st.metric("Total Datasets", len(datasets))

        if datasets:
            # Summary table
            summary_data = []
            for ds in datasets:
                qc = get_qc_dict(ds)
                summary_data.append(
                    {
                        "Name": ds.name,
                        "Omics Type": ds.omics_type,
                        "Ingestion Status": getattr(ds, "ingestion_status", "unknown"),
                        "Features": qc.get("num_features", ""),
                        "Missing Rate": (
                            f"{qc.get('missing_rate', ''):.3f}" if qc.get("missing_rate") is not None else ""
                        ),
                        "Source": getattr(ds, "external_ids", "") or getattr(ds, "file_paths", None) or "",
                        "Created": ds.created_at.strftime("%Y-%m-%d"),
                    }
                )
            df_summary = pd.DataFrame(summary_data)
            st.dataframe(df_summary, width='stretch', hide_index=True)

            # Export buttons
            col1, col2 = st.columns(2)
            with col1:
                csv_summary = df_summary.to_csv(index=False)
                st.download_button(
                    label="📥 Download Summary (CSV)",
                    data=csv_summary,
                    file_name="datasets_summary.csv",
                    mime="text/csv",
                )
            with col2:
                # Full dataset export
                full_data = []
                for ds in datasets:
                    full_data.append(
                        {
                            "ID": str(ds.id),
                            "Name": ds.name,
                            "Omics Type": ds.omics_type,
                            "Description": ds.description or "",
                            "Created": ds.created_at.strftime("%Y-%m-%d %H:%M:%S"),
                            "Updated": ds.updated_at.strftime("%Y-%m-%d %H:%M:%S"),
                            "File Paths": "; ".join(ds.file_paths) if ds.file_paths else "",
                            "Disease": ", ".join(ds.disease) if ds.disease else "",
                            "Notion ID": ds.notion_page_id or "",
                        }
                    )
                df_full = pd.DataFrame(full_data)
                csv_full = df_full.to_csv(index=False)
                st.download_button(
                    label="📥 Download Full Data (CSV)",
                    data=csv_full,
                    file_name="datasets_full.csv",
                    mime="text/csv",
                )

            st.markdown("---")
            st.subheader("Dataset Details")

            # Display datasets (access relationships while session is open)
            for dataset in datasets:
                with st.expander(f"**{dataset.name}** ({dataset.omics_type})"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{dataset.id}`")
                        st.write(f"**Omics Type:** {dataset.omics_type}")
                        if dataset.description:
                            st.write(f"**Description:** {dataset.description}")
                        qc = get_qc_dict(dataset)
                        ingestion_status = getattr(dataset, "ingestion_status", None)
                        if ingestion_status:  # --- Inline status badge ---
                            status_color = {"complete": "green", "in_progress": "orange", "failed": "red", "pending": "gray"}.get(ingestion_status, "gray")
                            st.markdown(
                                f"**Ingestion Status:** <span style='color:{status_color};font-weight:bold;'>{ingestion_status.upper()}</span>",
                                unsafe_allow_html=True,
                            )
                        if qc:
                            st.caption(
                                f"Data Summary: {qc.get('num_features', '?')} features, {qc.get('missing_rate', '?'):.1%} missing rate"
                            )
                    with col2:
                        st.write(f"**Created:** {dataset.created_at.strftime('%Y-%m-%d %H:%M')}")
                        st.write(f"**Updated:** {dataset.updated_at.strftime('%Y-%m-%d %H:%M')}")

                    if dataset.file_paths:
                        st.write("**File Paths:**")
                        for fp in dataset.file_paths:
                            st.code(fp, language=None)

                    if dataset.disease:
                        st.write(f"**Disease:** {', '.join(dataset.disease)}")

                    # Export as PowerPoint button
                    try:
                        pptx_data = generate_dataset_slides(dataset.id, db)
                        dataset_name_safe = (dataset.name or "dataset").replace(" ", "_").replace("/", "_")[:50]
                        st.download_button(
                            label="📊 Export as PowerPoint",
                            data=pptx_data,
                            file_name=f"{dataset_name_safe}.pptx",
                            mime_type="application/vnd.openxmlformats-officedocument.presentationml.presentation",
                            key=f"export_pptx_{dataset.id}",
                        )
                    except Exception as e:
                        st.error(f"Failed to generate PowerPoint: {e}")

                    # Related entities (Postgres relationships - access while session is open)
                    if dataset.programs:
                        st.write(f"**Linked Programs:** {len(dataset.programs)}")
                        program_names = [p.name for p in dataset.programs[:5]]
                        if program_names:
                            st.caption(
                                f"Programs: {', '.join(program_names)}"
                                + (f" ... and {len(dataset.programs) - 5} more" if len(dataset.programs) > 5 else "")
                            )

                    if dataset.experiments:
                        st.write(f"**Linked Experiments:** {len(dataset.experiments)}")
                        experiment_names = [e.name for e in dataset.experiments[:5]]
                        if experiment_names:
                            st.caption(
                                f"Experiments: {', '.join(experiment_names)}"
                                + (
                                    f" ... and {len(dataset.experiments) - 5} more"
                                    if len(dataset.experiments) > 5
                                    else ""
                                )
                            )

                    # Linked features (Postgres - access while session is open)
                    if dataset.features:
                        st.write(f"**Linked Features:** {len(dataset.features)}")
                        # Show feature types breakdown
                        feature_types = {}
                        for feature in dataset.features:
                            feature_types[feature.feature_type] = feature_types.get(feature.feature_type, 0) + 1
                        if feature_types:
                            type_breakdown = ", ".join([f"{ftype}: {count}" for ftype, count in feature_types.items()])
                            st.caption(f"Feature types: {type_breakdown}")

                        # Show sample features (first 10)
                        if len(dataset.features) > 0:
                            sample_features = dataset.features[:10]
                            feature_names = [f.name for f in sample_features]
                            st.caption(
                                f"Sample features: {', '.join(feature_names[:5])}"
                                + (f" ... and {len(dataset.features) - 5} more" if len(dataset.features) > 5 else "")
                            )

                    with st.expander("View Raw Metadata", expanded=False):
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
                    if dataset.experiments:
                        exp_names = ", ".join([e.name for e in dataset.experiments])
                        st.write(f"**Part of experiment(s):** {exp_names}")
        else:
            st.info("No datasets found matching your filters.")
