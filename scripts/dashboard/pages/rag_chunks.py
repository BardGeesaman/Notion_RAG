"""RAG Chunks page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy.orm import joinedload

from amprenta_rag.database.models import RAGChunk
from scripts.dashboard.db_session import db_session


def render_rag_chunks_page() -> None:
    """
    Render the RAG Chunks page with filtering, search, and detailed views.

    Features:
    - Filter by source type (Literature, Email, Dataset, etc.)
    - Search by chunk text, source name
    - Display source relationships
    - Export functionality (CSV)
    """
    st.header("📄 RAG Chunks")

    with db_session() as db:
        # Filters
        col1, col2 = st.columns(2)
        with col1:
            # Get distinct source types
            source_types = [row[0] for row in db.query(RAGChunk.source_type).distinct().all() if row[0]]
            source_type_filter = st.selectbox(
                "Filter by Source Type",
                ["All"] + source_types,
            )
        with col2:
            search_term = st.text_input("Search by chunk text or source name", "")

        # Query chunks with eager loading of relationships
        query = db.query(RAGChunk).options(
            joinedload(RAGChunk.literature),
            joinedload(RAGChunk.email),
        )

        if source_type_filter != "All":
            query = query.filter(RAGChunk.source_type == source_type_filter)

        if search_term:
            query = query.filter(
                (RAGChunk.chunk_text.ilike(f"%{search_term}%")) | (RAGChunk.source_name.ilike(f"%{search_term}%"))
            )

        chunks = query.order_by(RAGChunk.created_at.desc()).limit(100).all()

        st.metric("Chunks Found", len(chunks))
        if len(chunks) == 100:
            st.info("Showing first 100 chunks. Use filters to narrow down.")

        if chunks:
            # Summary table
            chunk_data = []
            for chunk in chunks:
                chunk_data.append(
                    {
                        "Chunk ID": chunk.chunk_id[:50] + "..." if len(chunk.chunk_id) > 50 else chunk.chunk_id,
                        "Source Type": chunk.source_type,
                        "Source Name": (
                            chunk.source_name[:50] + "..."
                            if chunk.source_name and len(chunk.source_name) > 50
                            else (chunk.source_name or "-")
                        ),
                        "Snippet": (
                            chunk.snippet[:60] + "..."
                            if chunk.snippet and len(chunk.snippet) > 60
                            else (chunk.snippet or "-")
                        ),
                        "Index": chunk.chunk_index,
                        "Created": chunk.created_at.strftime("%Y-%m-%d"),
                    }
                )
            df_chunks = pd.DataFrame(chunk_data)
            st.dataframe(df_chunks, use_container_width=True, hide_index=True)

            # Source type distribution
            if len(df_chunks) > 0:
                source_type_counts = df_chunks["Source Type"].value_counts()
                if len(source_type_counts) > 1:
                    st.subheader("Chunks by Source Type")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.bar_chart(source_type_counts)
                    with col2:
                        st.dataframe(
                            source_type_counts.reset_index().rename(
                                columns={"index": "Source Type", "Source Type": "Count"}
                            )
                        )

            # Export button
            csv_chunks = df_chunks.to_csv(index=False)
            st.download_button(
                label="📥 Download Chunks (CSV)",
                data=csv_chunks,
                file_name="rag_chunks.csv",
                mime="text/csv",
            )

            st.markdown("---")
            st.subheader("Chunk Details")

            for chunk in chunks:
                with st.expander(f"**Chunk {chunk.chunk_index}** from {chunk.source_name or chunk.source_type}"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**Chunk ID:** `{chunk.chunk_id}`")
                        st.write(f"**Source Type:** {chunk.source_type}")
                        if chunk.source_name:
                            st.write(f"**Source Name:** {chunk.source_name}")
                        if chunk.source_id:
                            st.write(f"**Source ID:** `{chunk.source_id}`")
                    with col2:
                        st.write(f"**Chunk Index:** {chunk.chunk_index}")
                        st.write(f"**Created:** {chunk.created_at.strftime('%Y-%m-%d %H:%M')}")
                        if chunk.zotero_item_key:
                            st.write(f"**Zotero Item Key:** `{chunk.zotero_item_key}`")

                    if chunk.snippet:
                        st.write("**Snippet:**")
                        st.text(chunk.snippet)

                    st.write("**Full Text:**")
                    st.text(chunk.chunk_text)

                    # Source relationships
                    if chunk.literature:
                        st.write(f"**Linked Literature:** {chunk.literature.title}")
                    if chunk.email:
                        st.write(f"**Linked Email:** {chunk.email.title}")
        else:
            st.info("No chunks found matching your filters.")
