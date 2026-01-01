"""Datasets page for the Streamlit dashboard."""

from __future__ import annotations

import json
import pandas as pd
import streamlit as st

from amprenta_rag.database.models import Dataset
from amprenta_rag.export.slide_generator import generate_dataset_slides
from amprenta_rag.notebooks import generate_dataset_notebook
from scripts.dashboard.db_session import db_session
from scripts.dashboard.components.comment_widget import render_comments_widget
from scripts.dashboard.components.annotation_panel import render_annotation_panel, render_annotation_indicator
from scripts.dashboard.utils.cache import fetch_datasets, clear_all_caches
from scripts.dashboard.utils.accessibility import (
    render_skip_link,
    add_heading_structure,
    ensure_minimum_contrast
)


def get_qc_dict(ds):
    return ds.qc_metrics if hasattr(ds, "qc_metrics") and ds.qc_metrics else {}


def render_datasets_page() -> None:
    """
    Render the Datasets page with filtering, search, and detailed views.

    Features:
    - Filter by omics type (cached data)
    - Search by name
    - Summary and detailed dataset views
    - Export functionality (CSV)
    - Cache refresh button
    """
    # Add accessibility features
    render_skip_link("main-datasets-content")
    ensure_minimum_contrast()
    
    # Add main content landmark and heading
    st.markdown(
        """
        <main id="main-datasets-content" role="main" aria-label="Dataset catalog and management">
        </main>
        """,
        unsafe_allow_html=True
    )
    
    add_heading_structure("📊 Datasets", level=1, id="datasets-title")
    
    # Add cache refresh button
    col_header, col_refresh = st.columns([10, 1])
    with col_header:
        st.markdown("")  # Spacing since we added our own heading
    with col_refresh:
        if st.button("🔄", help="Refresh data", key="refresh_datasets"):
            clear_all_caches()
            st.rerun()

    # Use cached data
    datasets_data = fetch_datasets(limit=200)
    
    # Extract distinct values for filters from cached data
    omics_types = list(set(d.get("omics_type") for d in datasets_data if d.get("omics_type")))
    omics_types.sort()
    
    status_opts = list(set(d.get("ingestion_status") for d in datasets_data if d.get("ingestion_status")))
    status_opts.sort()

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        omics_filter = st.selectbox(
            "Filter by Omics Type",
            ["All"] + omics_types,
        )
    with col2:
        qc_filter = st.selectbox("Ingestion Status", ["All"] + status_opts)
    with col3:
        search_term = st.text_input("Search by name", "")

    # Filter datasets (client-side filtering on cached data)
    filtered_datasets = datasets_data
    
    if omics_filter != "All":
        filtered_datasets = [d for d in filtered_datasets if d.get("omics_type") == omics_filter]
    
    if qc_filter != "All":
        filtered_datasets = [d for d in filtered_datasets if d.get("ingestion_status") == qc_filter]
    
    if search_term:
        filtered_datasets = [
            d for d in filtered_datasets 
            if search_term.lower() in d.get("name", "").lower()
        ]

    st.metric("Total Datasets", len(filtered_datasets))

    if filtered_datasets:
            # Summary table from cached data
            summary_data = []
            for ds_data in filtered_datasets:
                qc = ds_data.get("qc_metrics", {}) if isinstance(ds_data.get("qc_metrics"), dict) else {}
                summary_data.append(
                    {
                        "Name": ds_data.get("name", ""),
                        "Omics Type": ds_data.get("omics_type", ""),
                        "Ingestion Status": ds_data.get("ingestion_status", "unknown"),
                        "Features": qc.get("num_features", ""),
                        "Missing Rate": (
                            f"{qc.get('missing_rate', ''):.3f}" if qc.get("missing_rate") is not None else ""
                        ),
                        "Source": ds_data.get("external_ids", "") or ds_data.get("file_paths", "") or "",
                        "Created": ds_data.get("created_at", ""),
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
                # Full dataset export from cached data
                full_data = []
                for ds_data in filtered_datasets:
                    file_paths = ds_data.get("file_paths", [])
                    if isinstance(file_paths, list):
                        file_paths_str = "; ".join(file_paths)
                    else:
                        file_paths_str = str(file_paths) if file_paths else ""
                    
                    disease = ds_data.get("disease", [])
                    if isinstance(disease, list):
                        disease_str = ", ".join(disease)
                    else:
                        disease_str = str(disease) if disease else ""
                    
                    full_data.append(
                        {
                            "ID": str(ds_data.get("id", "")),
                            "Name": ds_data.get("name", ""),
                            "Omics Type": ds_data.get("omics_type", ""),
                            "Description": ds_data.get("description", ""),
                            "Created": ds_data.get("created_at", ""),
                            "Updated": ds_data.get("updated_at", ""),
                            "File Paths": file_paths_str,
                            "Disease": disease_str,
                            "Notion ID": ds_data.get("notion_page_id", ""),
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

            # Display datasets (use cached data for basic info, database for relationships)
            for ds_data in filtered_datasets:
                dataset_id = ds_data.get("id")
                # For complex operations like relationships, we still need database access
                with db_session() as db:
                    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                    if not dataset:
                        continue
                # Show annotation panel in sidebar if this dataset is selected
                if st.session_state.get("annotation_context", {}).get("entity_id") == str(dataset.id):
                    with st.sidebar:
                        render_annotation_panel(
                            entity_type="dataset",
                            entity_id=dataset.id,
                            position_type=st.session_state["annotation_context"].get("position_type"),
                            position_data=st.session_state["annotation_context"].get("position_data"),
                        )
                
                with st.expander(f"**{ds_data.get('name', 'Unknown')}** ({ds_data.get('omics_type', 'Unknown')})"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**ID:** `{ds_data.get('id')}`")
                        st.write(f"**Omics Type:** {ds_data.get('omics_type')}")
                        if ds_data.get("description"):
                            st.write(f"**Description:** {ds_data['description']}")
                        qc = ds_data.get("qc_metrics", {}) if isinstance(ds_data.get("qc_metrics"), dict) else {}
                        ingestion_status = ds_data.get("ingestion_status")
                        if ingestion_status:  # --- Inline status badge ---
                            status_color = {"complete": "green", "in_progress": "orange", "failed": "red", "pending": "gray"}.get(ingestion_status, "gray")
                            st.markdown(
                                f"**Ingestion Status:** <span style='color:{status_color};font-weight:bold;'>{ingestion_status.upper()}</span>",
                                unsafe_allow_html=True,
                            )
                        if qc:
                            missing_rate = qc.get('missing_rate')
                            missing_rate_str = f"{missing_rate:.1%}" if missing_rate is not None else "?"
                            st.caption(
                                f"Data Summary: {qc.get('num_features', '?')} features, {missing_rate_str} missing rate"
                            )
                    with col2:
                        st.write(f"**Created:** {ds_data.get('created_at', '')}")
                        st.write(f"**Updated:** {ds_data.get('updated_at', '')}")

                    file_paths = ds_data.get("file_paths", [])
                    if file_paths:
                        st.write("**File Paths:**")
                        if isinstance(file_paths, list):
                            for i, fp in enumerate(file_paths):
                                col_path, col_annotation = st.columns([10, 1])
                                with col_path:
                                    st.code(fp, language=None)
                            with col_annotation:
                                # Add annotation indicator for this file path
                                if render_annotation_indicator(
                                    entity_type="dataset",
                                    entity_id=dataset.id,
                                    position_type="field",
                                    position_data={"field": f"file_path_{i}"},
                                    label="📝"
                                ):
                                    pass

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

                    # Export as Jupyter Notebook (.ipynb)
                    try:
                        nb = generate_dataset_notebook(str(dataset.id))
                        dataset_name_safe = (dataset.name or "dataset").replace(" ", "_").replace("/", "_")[:50]
                        st.download_button(
                            label="📓 Export to Notebook",
                            data=json.dumps(nb),
                            file_name=f"{dataset_name_safe}.ipynb",
                            mime="application/json",
                            key=f"export_notebook_{dataset.id}",
                        )
                    except Exception as e:
                        st.error(f"Failed to generate notebook: {e}")

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

                        # Show sample features (first 10) with annotation indicators
                        if len(dataset.features) > 0:
                            st.write("**Sample Features:**")
                            sample_features = dataset.features[:5]  # Show fewer to make room for annotation buttons
                            
                            for feature in sample_features:
                                col_feature, col_annotation = st.columns([10, 1])
                                with col_feature:
                                    st.caption(f"• {feature.name} ({feature.feature_type})")
                                with col_annotation:
                                    # Add annotation indicator for this feature
                                    if render_annotation_indicator(
                                        entity_type="dataset",
                                        entity_id=dataset.id,
                                        position_type="column",
                                        position_data={"column": feature.name},
                                        label="📝"
                                    ):
                                        pass
                            
                            if len(dataset.features) > 5:
                                st.caption(f"... and {len(dataset.features) - 5} more features")

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
                    
                    # Comments section
                    st.markdown("---")
                    render_comments_widget("dataset", dataset.id)
    else:
        st.info("No datasets found matching your filters.")
