"""RAG Query page for the Streamlit dashboard."""

from __future__ import annotations

import streamlit as st

from amprenta_rag.query.rag_engine import query_rag


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

    # Query button
    if st.button("🔍 Search", type="primary"):
        if not query_text.strip():
            st.warning("Please enter a query.")
        else:
            with st.spinner("Searching..."):
                try:
                    # Build filters
                    filters = {}
                    if source_type_filter != "All":
                        filters["source_type"] = source_type_filter

                    # Query RAG
                    result = query_rag(
                        user_query=query_text,
                        top_k=top_k,
                        source_types=[source_type_filter] if source_type_filter != "All" else None,
                    )

                    if result and result.matches:
                        st.success(f"Found {len(result.matches)} results")

                        # Display answer if available
                        if result.answer:
                            st.subheader("Answer")
                            st.markdown(result.answer)
                            st.markdown("---")

                        # Display matches
                        st.subheader("Matches")
                        for i, match in enumerate(result.matches, 1):
                            with st.expander(f"**Result {i}** - Score: {match.score:.3f}"):
                                # Metadata
                                metadata = match.metadata or {}

                                col1, col2 = st.columns(2)
                                with col1:
                                    if metadata.get("source"):
                                        st.write(f"**Source:** {metadata['source']}")
                                    if metadata.get("source_type"):
                                        st.write(f"**Source Type:** {metadata['source_type']}")
                                    if metadata.get("title"):
                                        st.write(f"**Title:** {metadata['title']}")
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
