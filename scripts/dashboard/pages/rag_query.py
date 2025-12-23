"""RAG Query page for the Streamlit dashboard."""

from __future__ import annotations

import streamlit as st

from amprenta_rag.query.rag_engine import query_rag
from amprenta_rag.llm.model_registry import get_available_models
from amprenta_rag.llm.parallel_reasoning import parallel_query


def render_rag_query_page() -> None:
    """
    Render the RAG Query page for semantic search.

    Features:
    - Text query input
    - Filter by source type
    - Display query results with sources
    - Show chunk details
    """
    st.header("🔍 RAG Query")
    st.markdown("Query the RAG system using natural language to find relevant information across all data sources.")

    # Model selection
    available_models = get_available_models()
    model_options = {m["description"]: m["name"] for m in available_models}
    default_model = "gpt-4o"

    # Parallel mode checkbox
    use_parallel = st.checkbox("Use Multiple Models", value=False, help="Query multiple models in parallel and synthesize responses")

    selected_model_names = []  # Initialize for parallel mode
    selected_model_name = default_model  # Initialize for single mode

    if use_parallel:
        selected_models = st.multiselect(
            "Select models to use",
            options=list(model_options.keys()),
            default=list(model_options.keys())[:2] if len(model_options) >= 2 else list(model_options.keys()),
            help="Select which models to query in parallel"
        )
        selected_model_names = [model_options[m] for m in selected_models] if selected_models else []
    else:
        selected_model = st.selectbox(
            "Select model",
            options=list(model_options.keys()),
            index=list(model_options.values()).index(default_model) if default_model in model_options.values() else 0,
            help="Select the LLM model to use for querying"
        )
        selected_model_name = model_options[selected_model]

    # Query input
    query_text = st.text_area(
        "Enter your query",
        placeholder="e.g., What is the role of ceramide in ALS?",
        height=100,
    )

    # Filters
    col1, col2 = st.columns(2)
    with col1:
        top_k = st.slider("Number of results", min_value=1, max_value=50, value=10)
    with col2:
        source_type_filter = st.selectbox(
            "Filter by source type (optional)",
            ["All", "Literature", "Email", "Dataset", "Experiment"],
        )

    # Trust scoring option
    use_trust = st.checkbox(
        "Apply Trust Scoring",
        value=False,
        help="Weight results by source reliability (internal > external > general)"
    )

    use_graph_boost = st.checkbox(
        "Use Graph Boost",
        value=False,
        help="Boost RAG matches that are connected to query entities in the Evidence Graph",
    )

    # Query button
    if st.button("🔍 Search", type="primary"):
        if not query_text.strip():
            st.warning("Please enter a query.")
        elif use_parallel and not selected_model_names:
            st.warning("Please select at least one model for parallel mode.")
        else:
            with st.spinner("Searching..."):
                try:
                    if use_parallel:
                        # Parallel mode: query multiple models
                        if not selected_model_names:
                            st.error("No models selected for parallel query.")
                        else:
                            # First get RAG results (using default model for retrieval)
                            rag_result = query_rag(
                                user_query=query_text,
                                top_k=top_k,
                                source_types=[source_type_filter] if source_type_filter != "All" else None,
                                model=default_model,
                                use_trust_scoring=use_trust,
                                use_graph_boost=use_graph_boost,
                            )

                            # Then run parallel reasoning on the question
                            parallel_result = parallel_query(query_text, models=selected_model_names)

                            # Display RAG matches
                            if rag_result and rag_result.matches:
                                st.success(f"Found {len(rag_result.matches)} RAG results")

                                # Display trust summary if trust scoring enabled
                                if use_trust and rag_result.trust_summary:
                                    st.subheader("📊 Source Trust Analysis")
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        avg_trust = rag_result.trust_summary.get("average_trust", 0.0)
                                        st.metric("Average Trust", f"{avg_trust:.0%}")
                                    with col2:
                                        high_trust_count = rag_result.trust_summary.get("high_trust_count", 0)
                                        st.metric("High Trust Sources", high_trust_count)

                                    # Show trust level breakdown
                                    st.write("**Trust Level Distribution:**")
                                    levels = rag_result.trust_summary.get("levels", {})
                                    for level, count in levels.items():
                                        if count > 0:
                                            level_label = level.replace("_", " ").title()
                                            st.write(f"  - {level_label}: {count}")
                                    st.markdown("---")

                                # Display individual model responses
                                st.subheader("🤖 Individual Model Responses")
                                cols = st.columns(min(len(parallel_result["individual_responses"]), 3))

                                for idx, response in enumerate(parallel_result["individual_responses"]):
                                    col_idx = idx % len(cols)
                                    with cols[col_idx]:
                                        status_emoji = "✅" if response.get("success") else "❌"
                                        with st.expander(f"{status_emoji} {response['model']}", expanded=False):
                                            if response.get("success"):
                                                st.markdown(response.get("response", "No response"))
                                            else:
                                                st.error(f"Error: {response.get('error', 'Unknown error')}")

                                st.markdown("---")

                                # Display synthesized answer
                                st.subheader("✨ Synthesized Answer")
                                st.markdown(parallel_result["synthesis"])
                                st.markdown("---")

                                # Display RAG answer if available
                                if rag_result.answer:
                                    st.subheader("📚 RAG Context Answer")
                                    st.markdown(rag_result.answer)
                                    st.markdown("---")

                                # Display matches
                                st.subheader("Matches")
                                for i, match in enumerate(rag_result.matches, 1):
                                    # Get trust score if available
                                    metadata = match.metadata or {}
                                    trust_score = metadata.get("trust_score")
                                    if trust_score is not None:
                                        if trust_score >= 0.7:
                                            trust_badge = "🟢 High"
                                        elif trust_score >= 0.4:
                                            trust_badge = "🟡 Medium"
                                        else:
                                            trust_badge = "🔴 Low"
                                        score_display = f"**Result {i}** - Score: {match.score:.3f} {trust_badge}"
                                    else:
                                        score_display = f"**Result {i}** - Score: {match.score:.3f}"

                                    with st.expander(score_display):
                                        # Metadata
                                        col1, col2 = st.columns(2)
                                        with col1:
                                            if metadata.get("source"):
                                                source_text = f"**Source:** {metadata['source']}"
                                                if trust_badge and trust_score is not None:
                                                    source_text += f" {trust_badge}"
                                                st.write(source_text)
                                            if metadata.get("source_type"):
                                                st.write(f"**Source Type:** {metadata['source_type']}")
                                            if metadata.get("title"):
                                                st.write(f"**Title:** {metadata['title']}")
                                            if trust_score is not None:
                                                st.caption(f"Trust Score: {trust_score:.0%}")
                                        with col2:
                                            if metadata.get("chunk_id"):
                                                st.write(f"**Chunk ID:** `{metadata['chunk_id']}`")
                                            if metadata.get("chunk_index"):
                                                st.write(f"**Chunk Index:** {metadata['chunk_index']}")

                                        # Text content
                                        if match.text:
                                            st.write("**Content:**")
                                            st.text(match.text)

                                        # Additional metadata
                                        if metadata.get("diseases"):
                                            diseases = metadata["diseases"]
                                            if isinstance(diseases, list):
                                                st.caption(f"Diseases: {', '.join(diseases)}")
                                        if metadata.get("targets"):
                                            targets = metadata["targets"]
                                            if isinstance(targets, list):
                                                st.caption(f"Targets: {', '.join(targets)}")
                                        if metadata.get("features"):
                                            features = metadata["features"]
                                            if isinstance(features, list):
                                                st.caption(
                                                    f"Features: {', '.join(features[:10])}"
                                                    + (f" ... and {len(features) - 10} more" if len(features) > 10 else "")
                                                )
                            else:
                                st.info("No RAG results found. Try a different query or adjust filters.")
                    else:
                        # Single model mode: standard RAG query
                        # Build filters
                        filters = {}
                        if source_type_filter != "All":
                            filters["source_type"] = source_type_filter

                        # Query RAG
                        result = query_rag(
                            user_query=query_text,
                            top_k=top_k,
                            source_types=[source_type_filter] if source_type_filter != "All" else None,
                            model=selected_model_name,
                            use_trust_scoring=use_trust,
                            use_graph_boost=use_graph_boost,
                        )

                        if result and result.matches:
                            st.success(f"Found {len(result.matches)} results")

                            # Display trust summary if trust scoring enabled
                            if use_trust and result.trust_summary:
                                st.subheader("📊 Source Trust Analysis")
                                col1, col2 = st.columns(2)
                                with col1:
                                    avg_trust = result.trust_summary.get("average_trust", 0.0)
                                    st.metric("Average Trust", f"{avg_trust:.0%}")
                                with col2:
                                    high_trust_count = result.trust_summary.get("high_trust_count", 0)
                                    st.metric("High Trust Sources", high_trust_count)

                                # Show trust level breakdown
                                st.write("**Trust Level Distribution:**")
                                levels = result.trust_summary.get("levels", {})
                                for level, count in levels.items():
                                    if count > 0:
                                        level_label = level.replace("_", " ").title()
                                        st.write(f"  - {level_label}: {count}")
                                st.markdown("---")

                            # Display answer if available
                            if result.answer:
                                st.subheader("Answer")
                                st.markdown(result.answer)
                                st.markdown("---")

                            # Display matches
                            st.subheader("Matches")
                            for i, match in enumerate(result.matches, 1):
                                # Get trust score if available
                                metadata = match.metadata or {}
                                trust_score = metadata.get("trust_score")
                                if trust_score is not None:
                                    if trust_score >= 0.7:
                                        trust_badge = "🟢 High"
                                    elif trust_score >= 0.4:
                                        trust_badge = "🟡 Medium"
                                    else:
                                        trust_badge = "🔴 Low"
                                    score_display = f"**Result {i}** - Score: {match.score:.3f} {trust_badge}"
                                else:
                                    score_display = f"**Result {i}** - Score: {match.score:.3f}"

                                with st.expander(score_display):
                                    # Metadata
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        if metadata.get("source"):
                                            source_text = f"**Source:** {metadata['source']}"
                                            if trust_badge and trust_score is not None:
                                                source_text += f" {trust_badge}"
                                            st.write(source_text)
                                        if metadata.get("source_type"):
                                            st.write(f"**Source Type:** {metadata['source_type']}")
                                        if metadata.get("title"):
                                            st.write(f"**Title:** {metadata['title']}")
                                        if trust_score is not None:
                                            st.caption(f"Trust Score: {trust_score:.0%}")
                                    with col2:
                                        if metadata.get("chunk_id"):
                                            st.write(f"**Chunk ID:** `{metadata['chunk_id']}`")
                                        if metadata.get("chunk_index"):
                                            st.write(f"**Chunk Index:** {metadata['chunk_index']}")

                                    # Text content
                                    if match.text:
                                        st.write("**Content:**")
                                        st.text(match.text)

                                    # Additional metadata
                                    if metadata.get("diseases"):
                                        diseases = metadata["diseases"]
                                        if isinstance(diseases, list):
                                            st.caption(f"Diseases: {', '.join(diseases)}")
                                    if metadata.get("targets"):
                                        targets = metadata["targets"]
                                        if isinstance(targets, list):
                                            st.caption(f"Targets: {', '.join(targets)}")
                                    if metadata.get("features"):
                                        features = metadata["features"]
                                        if isinstance(features, list):
                                            st.caption(
                                                f"Features: {', '.join(features[:10])}"
                                                + (f" ... and {len(features) - 10} more" if len(features) > 10 else "")
                                            )
                        else:
                            st.info("No results found. Try a different query or adjust filters.")

                except Exception as e:
                    st.error(f"Error querying RAG: {str(e)}")
                    st.exception(e)

    # Example queries
    st.markdown("---")
    st.subheader("Example Queries")
    example_queries = [
        "What is the role of ceramide in neurodegeneration?",
        "Find datasets related to ALS",
        "What are the lipid signatures in CSF?",
        "Show me literature about sphingolipids",
    ]

    for example in example_queries:
        if st.button(f"💡 {example}", key=f"example_{example}"):
            st.session_state.query_text = example
            st.rerun()
