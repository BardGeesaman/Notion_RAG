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
    get_effect_size_preset,
)
from amprenta_rag.database.models import Dataset, Program, Signature
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
    tab1, tab2, tab3, tab4 = st.tabs(["Pathway Enrichment", "Program Signature Maps", "Signature Scoring", "Power Analysis"])

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

                                        # Download button
                                        csv_results = df_results.to_csv(index=False)
                                        st.download_button(
                                            label="📥 Download Results (CSV)",
                                            data=csv_results,
                                            file_name=f"pathway_enrichment_{dataset_id}.csv",
                                            mime="text/csv",
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

                                        csv_results = df_results.to_csv(index=False)
                                        st.download_button(
                                            label="📥 Download Results (CSV)",
                                            data=csv_results,
                                            file_name=f"pathway_enrichment_{signature_id}.csv",
                                            mime="text/csv",
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

                                csv_results = df_results.to_csv(index=False)
                                st.download_button(
                                    label="📥 Download Results (CSV)",
                                    data=csv_results,
                                    file_name="pathway_enrichment_custom.csv",
                                    mime="text/csv",
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

                                        # Download results
                                        csv = df_matches.to_csv(index=False)
                                        st.download_button(
                                            label="📥 Download Results (CSV)",
                                            data=csv,
                                            file_name=f"signature_scores_{signature_id}_{dataset_id}.csv",
                                            mime="text/csv",
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
