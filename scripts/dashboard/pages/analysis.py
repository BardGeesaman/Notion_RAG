"""Analysis Tools page for the Streamlit dashboard."""

from __future__ import annotations

from uuid import UUID

import pandas as pd
import streamlit as st
import plotly.graph_objects as go

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.analysis.program_maps.reporting import generate_program_map_report, generate_program_signature_map
from amprenta_rag.analysis.power_analysis import (
    calculate_sample_size,
    calculate_power,
)
from amprenta_rag.analysis.pattern_detection import (
    find_recurring_features,
    get_cross_dataset_summary,
)
from amprenta_rag.analysis.confounder_detection import (
    detect_confounders,
    get_confounder_report,
)
from amprenta_rag.analysis.design_engine import (
    recommend_design,
    get_design_requirements,
    validate_design,
)
from amprenta_rag.database.models import Dataset, Program, Signature
from amprenta_rag.utils.data_export import export_to_csv, export_to_json, export_to_excel
from scripts.dashboard.db_session import db_session


def render_analysis_page() -> None:
    """
    Render the Analysis Tools page for pathway enrichment, program maps, and signature scoring.

    Features:
    - Pathway enrichment analysis
    - Program signature maps
    - Signature scoring against datasets
    """
    st.header("🔬 Analysis Tools")
    st.markdown("Perform pathway enrichment, program signature mapping, and signature scoring.")

    # Tab selection
    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(["Pathway Enrichment", "Program Signature Maps", "Signature Scoring", "Power Analysis", "Cross-Dataset Patterns", "Confounder Check", "Design Assistant"])

    # Tab 1: Pathway Enrichment
    with tab1:
        st.subheader("Pathway Enrichment Analysis")
        st.markdown("Identify pathways significantly enriched in your datasets or signatures.")

        analysis_type = st.radio("Select Analysis Type", ["Dataset", "Signature", "Custom Features"], horizontal=True)

        if analysis_type == "Dataset":
            with db_session() as db:
                datasets = db.query(Dataset).order_by(Dataset.name).all()

                if datasets:
                    dataset_options = {f"{d.name} ({d.id})": str(d.id) for d in datasets}
                    selected_dataset = st.selectbox(
                        "Select a dataset",
                        list(dataset_options.keys()),
                    )
                    dataset_id = UUID(dataset_options[selected_dataset])

                    # Get dataset features
                    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                    if dataset and dataset.features:
                        feature_names = {f.name for f in dataset.features}
                        feature_types = {f.feature_type for f in dataset.features}

                        st.info(f"Found {len(feature_names)} features in dataset")

                        # Parameters
                        col1, col2 = st.columns(2)
                        with col1:
                            p_value_threshold = st.slider(
                                "P-value Threshold",
                                min_value=0.001,
                                max_value=0.1,
                                value=0.05,
                                step=0.001,
                                format="%.3f",
                            )
                            top_pathways = st.number_input(
                                "Top Pathways to Display", min_value=5, max_value=50, value=10
                            )
                        with col2:
                            pathway_sources = st.multiselect(
                                "Pathway Sources", ["KEGG", "Reactome"], default=["KEGG", "Reactome"]
                            )

                        if st.button("🔬 Run Enrichment Analysis", type="primary"):
                            with st.spinner("Running pathway enrichment analysis... This may take a few minutes."):
                                try:
                                    results = perform_pathway_enrichment(
                                        input_features=feature_names,
                                        input_feature_types=feature_types,
                                        pathway_sources=pathway_sources if pathway_sources else None,
                                        p_value_threshold=p_value_threshold,
                                    )

                                    if results:
                                        st.success(f"✅ Found {len(results)} significantly enriched pathways!")

                                        # Display results table
                                        results_data = []
                                        for r in results[:top_pathways]:
                                            results_data.append(
                                                {
                                                    "Pathway": r.pathway.name,
                                                    "Source": r.pathway.source,
                                                    "ID": r.pathway.id,
                                                    "Input Features": r.input_features_in_pathway,
                                                    "Pathway Size": r.pathway_size,
                                                    "P-value": f"{r.p_value:.4f}",
                                                    "Adjusted P-value": f"{r.adjusted_p_value:.4f}",
                                                    "Odds Ratio": f"{r.odds_ratio:.2f}" if r.odds_ratio else "N/A",
                                                }
                                            )

                                        df_results = pd.DataFrame(results_data)
                                        st.dataframe(df_results, width='stretch', hide_index=True)

                                        # Export section
                                        st.markdown("### Export")
                                        col1, col2 = st.columns([1, 3])
                                        with col1:
                                            export_format = st.selectbox("Format", ["csv", "json", "excel"], key="pathway_export_format_ds")
                                        with col2:
                                            mime_types = {"csv": "text/csv", "json": "application/json", "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}
                                            file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                                            if export_format == "csv":
                                                export_data = export_to_csv(df_results)
                                            elif export_format == "json":
                                                export_data = export_to_json(df_results)
                                            else:
                                                export_data = export_to_excel(df_results)

                                            st.download_button(
                                                label=f"📥 Download Results ({export_format.upper()})",
                                                data=export_data,
                                                file_name=f"pathway_enrichment_{dataset_id}.{file_extensions[export_format]}",
                                                mime=mime_types[export_format],
                                            )
                                    else:
                                        st.warning(
                                            "No significantly enriched pathways found. Try adjusting the p-value threshold."
                                        )

                                except Exception as e:
                                    st.error(f"❌ Enrichment analysis failed: {str(e)}")
                                    st.exception(e)
                    else:
                        st.warning("Dataset has no features. Cannot perform enrichment analysis.")
                else:
                    st.info("No datasets available. Ingest a dataset first.")

        elif analysis_type == "Signature":
            with db_session() as db:
                signatures = db.query(Signature).order_by(Signature.name).all()

                if signatures:
                    signature_options = {f"{s.name} ({s.id})": str(s.id) for s in signatures}
                    selected_signature = st.selectbox(
                        "Select a signature",
                        list(signature_options.keys()),
                    )
                    signature_id = UUID(signature_options[selected_signature])

                    # Get signature features
                    signature = db.query(Signature).filter(Signature.id == signature_id).first()
                    if signature and signature.features:
                        feature_names = {f.name for f in signature.features}
                        feature_types = {f.feature_type for f in signature.features}

                        st.info(f"Found {len(feature_names)} features in signature")

                        # Parameters (same as dataset)
                        col1, col2 = st.columns(2)
                        with col1:
                            p_value_threshold = st.slider(
                                "P-value Threshold",
                                min_value=0.001,
                                max_value=0.1,
                                value=0.05,
                                step=0.001,
                                format="%.3f",
                                key="sig_pval",
                            )
                            top_pathways = st.number_input(
                                "Top Pathways to Display", min_value=5, max_value=50, value=10, key="sig_top"
                            )
                        with col2:
                            pathway_sources = st.multiselect(
                                "Pathway Sources", ["KEGG", "Reactome"], default=["KEGG", "Reactome"], key="sig_sources"
                            )

                        if st.button("🔬 Run Enrichment Analysis", type="primary", key="sig_btn"):
                            with st.spinner("Running pathway enrichment analysis... This may take a few minutes."):
                                try:
                                    results = perform_pathway_enrichment(
                                        input_features=feature_names,
                                        input_feature_types=feature_types,
                                        pathway_sources=pathway_sources if pathway_sources else None,
                                        p_value_threshold=p_value_threshold,
                                    )

                                    if results:
                                        st.success(f"✅ Found {len(results)} significantly enriched pathways!")

                                        # Display results (same format as dataset)
                                        results_data = []
                                        for r in results[:top_pathways]:
                                            results_data.append(
                                                {
                                                    "Pathway": r.pathway.name,
                                                    "Source": r.pathway.source,
                                                    "ID": r.pathway.id,
                                                    "Input Features": r.input_features_in_pathway,
                                                    "Pathway Size": r.pathway_size,
                                                    "P-value": f"{r.p_value:.4f}",
                                                    "Adjusted P-value": f"{r.adjusted_p_value:.4f}",
                                                    "Odds Ratio": f"{r.odds_ratio:.2f}" if r.odds_ratio else "N/A",
                                                }
                                            )

                                        df_results = pd.DataFrame(results_data)
                                        st.dataframe(df_results, width='stretch', hide_index=True)

                                        # Export section
                                        st.markdown("### Export")
                                        col1, col2 = st.columns([1, 3])
                                        with col1:
                                            export_format = st.selectbox("Format", ["csv", "json", "excel"], key="pathway_export_format_sig")
                                        with col2:
                                            mime_types = {"csv": "text/csv", "json": "application/json", "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}
                                            file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                                            if export_format == "csv":
                                                export_data = export_to_csv(df_results)
                                            elif export_format == "json":
                                                export_data = export_to_json(df_results)
                                            else:
                                                export_data = export_to_excel(df_results)

                                            st.download_button(
                                                label=f"📥 Download Results ({export_format.upper()})",
                                                data=export_data,
                                                file_name=f"pathway_enrichment_{signature_id}.{file_extensions[export_format]}",
                                                mime=mime_types[export_format],
                                            )
                                    else:
                                        st.warning("No significantly enriched pathways found.")

                                except Exception as e:
                                    st.error(f"❌ Enrichment analysis failed: {str(e)}")
                                    st.exception(e)
                    else:
                        st.warning("Signature has no features. Cannot perform enrichment analysis.")
                else:
                    st.info("No signatures available. Ingest a signature first.")

        else:  # Custom Features
            st.markdown("### Enter Custom Features")
            feature_input = st.text_area(
                "Enter feature names (one per line or comma-separated)",
                height=150,
                help="Enter feature identifiers like gene names, metabolite names, etc.",
            )

            feature_type = st.selectbox(
                "Feature Type",
                ["gene", "protein", "metabolite", "lipid"],
                help="Select the type of features you're entering",
            )

            if feature_input:
                # Parse features
                features = set()
                for line in feature_input.split("\n"):
                    for item in line.split(","):
                        item = item.strip()
                        if item:
                            features.add(item)

                st.info(f"Parsed {len(features)} features")

                # Parameters
                col1, col2 = st.columns(2)
                with col1:
                    p_value_threshold = st.slider(
                        "P-value Threshold",
                        min_value=0.001,
                        max_value=0.1,
                        value=0.05,
                        step=0.001,
                        format="%.3f",
                        key="custom_pval",
                    )
                    top_pathways = st.number_input(
                        "Top Pathways to Display", min_value=5, max_value=50, value=10, key="custom_top"
                    )
                with col2:
                    pathway_sources = st.multiselect(
                        "Pathway Sources", ["KEGG", "Reactome"], default=["KEGG", "Reactome"], key="custom_sources"
                    )

                if st.button("🔬 Run Enrichment Analysis", type="primary", key="custom_btn"):
                    with st.spinner("Running pathway enrichment analysis... This may take a few minutes."):
                        try:
                            results = perform_pathway_enrichment(
                                input_features=features,
                                input_feature_types={feature_type},
                                pathway_sources=pathway_sources if pathway_sources else None,
                                p_value_threshold=p_value_threshold,
                            )

                            if results:
                                st.success(f"✅ Found {len(results)} significantly enriched pathways!")

                                # Display results
                                results_data = []
                                for r in results[:top_pathways]:
                                    results_data.append(
                                        {
                                            "Pathway": r.pathway.name,
                                            "Source": r.pathway.source,
                                            "ID": r.pathway.id,
                                            "Input Features": r.input_features_in_pathway,
                                            "Pathway Size": r.pathway_size,
                                            "P-value": f"{r.p_value:.4f}",
                                            "Adjusted P-value": f"{r.adjusted_p_value:.4f}",
                                            "Odds Ratio": f"{r.odds_ratio:.2f}" if r.odds_ratio else "N/A",
                                        }
                                    )

                                df_results = pd.DataFrame(results_data)
                                st.dataframe(df_results, width='stretch', hide_index=True)

                                # Export section
                                st.markdown("### Export")
                                col1, col2 = st.columns([1, 3])
                                with col1:
                                    export_format = st.selectbox("Format", ["csv", "json", "excel"], key="pathway_export_format_custom")
                                with col2:
                                    mime_types = {"csv": "text/csv", "json": "application/json", "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}
                                    file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                                    if export_format == "csv":
                                        export_data = export_to_csv(df_results)
                                    elif export_format == "json":
                                        export_data = export_to_json(df_results)
                                    else:
                                        export_data = export_to_excel(df_results)

                                    st.download_button(
                                        label=f"📥 Download Results ({export_format.upper()})",
                                        data=export_data,
                                        file_name=f"pathway_enrichment_custom.{file_extensions[export_format]}",
                                        mime=mime_types[export_format],
                                    )
                            else:
                                st.warning("No significantly enriched pathways found.")

                        except Exception as e:
                            st.error(f"❌ Enrichment analysis failed: {str(e)}")
                            st.exception(e)

    # Tab 2: Program Signature Maps
    with tab2:
        st.subheader("Program Signature Maps")
        st.markdown("Generate program × signature scoring matrices and coverage analysis.")

        with db_session() as db:
            programs = db.query(Program).order_by(Program.name).all()

            if programs:
                program_options = {f"{p.name} ({p.id})": str(p.id) for p in programs}
                selected_program = st.selectbox(
                    "Select a program",
                    list(program_options.keys()),
                )
                program_id = UUID(program_options[selected_program])

                col1, col2 = st.columns(2)
                with col1:
                    top_n = st.number_input(
                        "Top N Signatures",
                        min_value=5,
                        max_value=50,
                        value=10,
                        help="Number of top signatures to include in the map",
                    )

                if st.button("🗺️ Generate Program Map", type="primary"):
                    with st.spinner("Generating program signature map... This may take a few minutes."):
                        try:
                            # Generate map
                            program_map = generate_program_signature_map(
                                program_page_id=str(program_id),  # Note: may need UUID conversion
                                top_n=top_n,
                            )

                            # Generate report
                            report = generate_program_map_report(program_map)

                            st.success("✅ Program signature map generated!")
                            st.markdown("---")
                            st.markdown("### Map Report")
                            st.markdown(report)

                            # Download report
                            st.download_button(
                                label="📥 Download Report (Markdown)",
                                data=report,
                                file_name=f"program_map_{program_id}.md",
                                mime="text/markdown",
                            )

                        except Exception as e:
                            st.error(f"❌ Program map generation failed: {str(e)}")
                            st.exception(e)
            else:
                st.info("No programs available. Create a program first.")

    # Tab 3: Signature Scoring
    with tab3:
        st.subheader("Signature Scoring")
        st.markdown("Score a signature against a dataset to measure match quality.")

        with db_session() as db:
            signatures = db.query(Signature).order_by(Signature.name).all()
            datasets = db.query(Dataset).order_by(Dataset.name).all()

            if signatures and datasets:
                col1, col2 = st.columns(2)
                with col1:
                    signature_options = {f"{s.name} ({s.id})": str(s.id) for s in signatures}
                    selected_signature = st.selectbox(
                        "Select a signature",
                        list(signature_options.keys()),
                    )
                    signature_id = UUID(signature_options[selected_signature])

                with col2:
                    dataset_options = {f"{d.name} ({d.id})": str(d.id) for d in datasets}
                    selected_dataset = st.selectbox(
                        "Select a dataset",
                        list(dataset_options.keys()),
                    )
                    dataset_id = UUID(dataset_options[selected_dataset])

                if st.button("📊 Score Signature", type="primary"):
                    with st.spinner("Scoring signature against dataset..."):
                        try:
                            from amprenta_rag.database.models import Dataset as DatasetModel
                            from amprenta_rag.ingestion.postgres_signature_matching import (
                                find_matching_signatures_for_postgres_dataset,
                            )

                            with db_session() as score_db:
                                # Get dataset
                                dataset = score_db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()

                                if not dataset:
                                    st.error("Dataset not found.")
                                else:
                                    # Score signature
                                    matches = find_matching_signatures_for_postgres_dataset(
                                        dataset_id=dataset_id,
                                        overlap_threshold=0.3,
                                        omics_type=dataset.omics_type,
                                    )

                                    if matches:
                                        st.success(f"Found {len(matches)} matching signature(s)!")

                                        # Display results
                                        st.subheader("Match Results")
                                        match_data = []
                                        for match in matches:
                                            # SignatureMatchResult has signature_name, score, overlap_fraction, etc.
                                            signature_name = getattr(match, "signature_name", "Unknown")
                                            if signature_name == "Unknown":
                                                # Try to get from database using signature_page_id (UUID)
                                                try:
                                                    from amprenta_rag.database.models import Signature as SignatureModel

                                                    sig_model = (
                                                        score_db.query(SignatureModel)
                                                        .filter(SignatureModel.id == UUID(match.signature_page_id))
                                                        .first()
                                                    )
                                                    if sig_model:
                                                        signature_name = sig_model.name
                                                except Exception:
                                                    pass

                                            # Get component counts from score_result if available
                                            matched_count = len(getattr(match, "matched_components", []))
                                            total_components = matched_count + len(
                                                getattr(match, "missing_components", [])
                                            )

                                            match_data.append(
                                                {
                                                    "Signature": signature_name,
                                                    "Overlap Fraction": f"{match.overlap_fraction:.3f}",
                                                    "Score": f"{match.score:.3f}",
                                                    "Matched Components": matched_count,
                                                    "Total Components": (
                                                        total_components if total_components > 0 else "N/A"
                                                    ),
                                                }
                                            )

                                        df_matches = pd.DataFrame(match_data)
                                        st.dataframe(df_matches, width='stretch', hide_index=True)

                                        # Export section
                                        st.markdown("### Export")
                                        col1, col2 = st.columns([1, 3])
                                        with col1:
                                            export_format = st.selectbox("Format", ["csv", "json", "excel"], key="signature_score_export_format")
                                        with col2:
                                            mime_types = {"csv": "text/csv", "json": "application/json", "excel": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"}
                                            file_extensions = {"csv": "csv", "json": "json", "excel": "xlsx"}

                                            if export_format == "csv":
                                                export_data = export_to_csv(df_matches)
                                            elif export_format == "json":
                                                export_data = export_to_json(df_matches)
                                            else:
                                                export_data = export_to_excel(df_matches)

                                            st.download_button(
                                                label=f"📥 Download Results ({export_format.upper()})",
                                                data=export_data,
                                                file_name=f"signature_scores_{signature_id}_{dataset_id}.{file_extensions[export_format]}",
                                                mime=mime_types[export_format],
                                            )
                                    else:
                                        st.info("No matching signatures found (overlap threshold: 0.3)")
                        except Exception as e:
                            st.error(f"❌ Error scoring signature: {str(e)}")
                            st.exception(e)
            else:
                if not signatures:
                    st.warning("No signatures available. Ingest a signature first.")
                if not datasets:
                    st.warning("No datasets available. Ingest a dataset first.")

    # Tab 4: Power Analysis
    with tab4:
        st.subheader("Statistical Power Analysis")
        st.markdown("Calculate required sample sizes and statistical power for your experiments.")

        col1, col2 = st.columns(2)

        with col1:
            test_type = st.selectbox(
                "Test Type",
                ["t-test", "anova", "correlation", "chi-square"],
                help="Select the statistical test you plan to use",
            )

            st.markdown("**Effect Size**")
            effect_size_col1, effect_size_col2, effect_size_col3 = st.columns(3)
            with effect_size_col1:
                if st.button("Small (0.2)", use_container_width=True):
                    st.session_state["effect_size"] = 0.2
            with effect_size_col2:
                if st.button("Medium (0.5)", use_container_width=True):
                    st.session_state["effect_size"] = 0.5
            with effect_size_col3:
                if st.button("Large (0.8)", use_container_width=True):
                    st.session_state["effect_size"] = 0.8

            effect_size = st.number_input(
                "Effect Size",
                min_value=0.01,
                max_value=2.0,
                value=st.session_state.get("effect_size", 0.5),
                step=0.01,
                help="Cohen's d for t-test, f for ANOVA, r for correlation",
            )

        with col2:
            alpha = st.slider(
                "Alpha (Significance Level)",
                min_value=0.01,
                max_value=0.10,
                value=0.05,
                step=0.01,
                help="Type I error rate",
            )

            power = st.slider(
                "Power",
                min_value=0.70,
                max_value=0.95,
                value=0.80,
                step=0.05,
                help="Statistical power (1 - Type II error rate)",
            )

        if st.button("Calculate Sample Size", type="primary"):
            try:
                required_n = calculate_sample_size(
                    effect_size=effect_size,
                    alpha=alpha,
                    power=power,
                    test_type=test_type,
                )

                st.success(f"✅ Required sample size: **{required_n}** {'per group' if test_type == 't-test' else 'total'}")

                # Power vs Sample Size chart
                st.markdown("### Power vs Sample Size")
                sample_sizes = list(range(max(5, required_n - 50), required_n + 100, 5))
                powers = [
                    calculate_power(n=n, effect_size=effect_size, alpha=alpha, test_type=test_type)
                    for n in sample_sizes
                ]

                fig = go.Figure()
                fig.add_trace(
                    go.Scatter(
                        x=sample_sizes,
                        y=powers,
                        mode="lines",
                        name="Power",
                        line=dict(color="blue", width=2),
                    )
                )
                fig.add_hline(
                    y=power,
                    line_dash="dash",
                    line_color="red",
                    annotation_text=f"Target Power ({power})",
                )
                fig.add_vline(
                    x=required_n,
                    line_dash="dash",
                    line_color="green",
                    annotation_text=f"Required N ({required_n})",
                )
                fig.update_layout(
                    title="Power vs Sample Size",
                    xaxis_title="Sample Size",
                    yaxis_title="Power",
                    height=400,
                )
                st.plotly_chart(fig, use_container_width=True)

            except Exception as e:
                st.error(f"❌ Error calculating sample size: {str(e)}")
                st.exception(e)

    # Tab 5: Cross-Dataset Patterns
    with tab5:
        st.subheader("Cross-Dataset Pattern Detection")
        st.markdown("Find recurring features and patterns across multiple datasets.")

        with db_session() as db:
            datasets = db.query(Dataset).order_by(Dataset.name).all()

            if not datasets:
                st.info("No datasets available. Ingest datasets first.")
            else:
                dataset_options = {f"{d.name} ({d.id})": str(d.id) for d in datasets}
                selected_datasets = st.multiselect(
                    "Select datasets for pattern analysis",
                    list(dataset_options.keys()),
                    help="Select 2 or more datasets to find recurring features",
                )

                if len(selected_datasets) >= 2:
                    dataset_ids = [dataset_options[d] for d in selected_datasets]

                    col1, col2 = st.columns(2)
                    with col1:
                        min_occurrence = st.number_input(
                            "Minimum Occurrence",
                            min_value=2,
                            max_value=len(selected_datasets),
                            value=2,
                            help="Minimum number of datasets a feature must appear in",
                        )

                    with col2:
                        if st.button("Find Patterns", type="primary"):
                            with st.spinner("Analyzing patterns across datasets..."):
                                try:
                                    # Find recurring features
                                    recurring = find_recurring_features(
                                        dataset_ids=dataset_ids,
                                        db=db,
                                        min_occurrence=min_occurrence,
                                    )

                                    if recurring:
                                        st.markdown("### Recurring Features")
                                        recurring_df = pd.DataFrame([
                                            {
                                                "Feature": r["feature_name"],
                                                "Occurrence Count": r["occurrence_count"],
                                                "Datasets": ", ".join([d[:8] + "..." for d in r["datasets"]]),
                                            }
                                            for r in recurring
                                        ])
                                        st.dataframe(recurring_df, hide_index=True, use_container_width=True)

                                        # Summary statistics
                                        st.markdown("### Summary Statistics")
                                        summary = get_cross_dataset_summary(dataset_ids, db)

                                        col_a, col_b = st.columns(2)
                                        with col_a:
                                            st.markdown("**Total Features per Dataset:**")
                                            for did, count in summary["total_features_per_dataset"].items():
                                                dataset_name = next(
                                                    (d.name for d in datasets if str(d.id) == did),
                                                    did[:8] + "..."
                                                )
                                                st.metric(dataset_name, count)

                                        with col_b:
                                            st.markdown("**Unique Features per Dataset:**")
                                            for did, count in summary["unique_features_per_dataset"].items():
                                                dataset_name = next(
                                                    (d.name for d in datasets if str(d.id) == did),
                                                    did[:8] + "..."
                                                )
                                                st.metric(f"{dataset_name} (unique)", count)

                                        if summary["overlap_counts"]:
                                            st.markdown("**Pairwise Overlaps:**")
                                            overlap_df = pd.DataFrame([
                                                {"Dataset Pair": pair, "Overlap Count": count}
                                                for pair, count in summary["overlap_counts"].items()
                                            ])
                                            st.dataframe(overlap_df, hide_index=True, use_container_width=True)
                                    else:
                                        st.info(f"No features found recurring in {min_occurrence} or more datasets.")

                                except Exception as e:
                                    st.error(f"❌ Error analyzing patterns: {str(e)}")
                                    st.exception(e)
                else:
                    st.info("Select at least 2 datasets to perform pattern analysis.")

    # Tab 6: Confounder Check
    with tab6:
        st.subheader("Confounder Detection")
        st.markdown("Upload sample metadata and detect potential confounders in experimental design.")

        uploaded_file = st.file_uploader(
            "Upload Sample Metadata CSV",
            type=["csv"],
            help="CSV file with sample metadata (one row per sample)",
        )

        if uploaded_file is not None:
            try:
                metadata_df = pd.read_csv(uploaded_file)
                st.success(f"Loaded {len(metadata_df)} samples with {len(metadata_df.columns)} columns")

                # Display preview
                with st.expander("Preview Data"):
                    st.dataframe(metadata_df.head(10), use_container_width=True)

                # Select group column
                group_column = st.selectbox(
                    "Select Group Column",
                    metadata_df.columns.tolist(),
                    help="Column containing group labels (e.g., treatment vs control)",
                )

                if st.button("Detect Confounders", type="primary"):
                    with st.spinner("Analyzing potential confounders..."):
                        try:
                            # Detect confounders
                            results = detect_confounders(metadata_df, group_column)

                            if results:
                                # Display results table
                                st.markdown("### Confounder Detection Results")
                                results_df = pd.DataFrame(results)
                                results_df["p_value"] = results_df["p_value"].apply(lambda x: f"{x:.4f}" if x is not None else "N/A")
                                results_df["is_confounder"] = results_df["is_confounder"].apply(lambda x: "⚠️ Yes" if x else "✓ No")
                                results_df.columns = ["Column", "Test", "P-Value", "Confounder"]
                                st.dataframe(results_df, hide_index=True, use_container_width=True)

                                # Get report and show warnings
                                report = get_confounder_report(metadata_df, group_column)

                                st.markdown("### Summary")
                                st.info(report["summary"])

                                if report["warnings"]:
                                    st.markdown("### Warnings")
                                    for warning in report["warnings"]:
                                        st.warning(warning)

                                # Highlight significant confounders
                                confounders = [r for r in results if r["is_confounder"]]
                                if confounders:
                                    st.markdown("### ⚠️ Significant Confounders Detected")
                                    confounder_df = pd.DataFrame([
                                        {
                                            "Column": c["column"],
                                            "Test": c["test"],
                                            "P-Value": f"{c['p_value']:.4f}",
                                        }
                                        for c in confounders
                                    ])
                                    st.dataframe(confounder_df, hide_index=True, use_container_width=True)
                                    st.markdown(
                                        "**Recommendation:** Consider stratifying or adjusting for these variables "
                                        "in your analysis to avoid confounding effects."
                                    )
                            else:
                                st.info("No variables tested. Check that your metadata file has valid data.")

                        except Exception as e:
                            st.error(f"❌ Error detecting confounders: {str(e)}")
                            st.exception(e)

            except Exception as e:
                st.error(f"❌ Error reading CSV file: {str(e)}")
                st.exception(e)
        else:
            st.info("Upload a CSV file with sample metadata to begin confounder detection.")

    # Tab 7: Design Assistant
    with tab7:
        st.subheader("Experimental Design Assistant")
        st.markdown("Get AI-powered recommendations for your experimental design based on your research question.")

        research_question = st.text_area(
            "Research Question",
            height=100,
            placeholder="e.g., Does treatment X reduce disease severity compared to placebo?",
            help="Describe what you want to investigate",
        )

        col1, col2 = st.columns(2)
        with col1:
            sample_count = st.number_input(
                "Available Samples",
                min_value=1,
                value=20,
                help="Total number of samples you have or plan to collect",
            )

        with col2:
            variables_input = st.text_input(
                "Variables/Factors",
                placeholder="e.g., treatment, age, gender",
                help="Comma-separated list of variables or factors",
            )

        if st.button("Get Recommendations", type="primary"):
            if not research_question:
                st.error("Please enter a research question")
            else:
                with st.spinner("Analyzing research question and generating recommendations..."):
                    try:
                        # Parse variables
                        variables = [v.strip() for v in variables_input.split(",")] if variables_input else []

                        # Get recommendation
                        recommendation = recommend_design(
                            research_question=research_question,
                            sample_count=sample_count,
                            variables=variables,
                        )

                        design_type = recommendation.get("design_type", "observational")

                        # Display recommendation
                        st.markdown("### Recommended Design")
                        st.success(f"**Design Type:** {design_type.replace('_', ' ').title()}")

                        st.markdown("**Rationale:**")
                        st.write(recommendation.get("rationale", "No rationale provided"))

                        # Get requirements
                        requirements = get_design_requirements(design_type)

                        st.markdown("### Design Requirements")
                        st.info(f"**Description:** {requirements['description']}")

                        st.markdown("**Minimum Requirements:**")
                        for req in requirements["requirements"]:
                            st.markdown(f"- {req}")

                        st.metric("Minimum Samples", requirements["min_samples"])
                        st.metric("Minimum Groups", requirements["min_groups"])

                        # Validate design
                        estimated_groups = len(variables) + 1 if variables else 2  # Rough estimate
                        warnings = validate_design(design_type, sample_count, estimated_groups)

                        if warnings:
                            st.markdown("### ⚠️ Design Validation Warnings")
                            for warning in warnings:
                                st.warning(warning)
                        else:
                            st.success("✓ Design meets minimum requirements")

                        # Show considerations
                        considerations = recommendation.get("considerations", [])
                        if considerations:
                            st.markdown("### Additional Considerations")
                            for consideration in considerations:
                                st.markdown(f"- {consideration}")

                    except Exception as e:
                        st.error(f"❌ Error generating recommendations: {str(e)}")
                        st.exception(e)
