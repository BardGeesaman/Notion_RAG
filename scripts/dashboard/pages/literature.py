"""Literature page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy.orm import joinedload

from amprenta_rag.database.models import Literature
from scripts.dashboard.db_session import db_session


def render_literature_page() -> None:
    """
    Render the Literature page with filtering, search, and detailed views.

    Features:
    - Filter by source type (journalArticle, book, etc.)
    - Search by title, author, DOI
    - Display embedding status
    - Show linked RAG chunks
    - Export functionality (CSV)
    """
    st.header("📚 Literature")

    with db_session() as db:
        # Filters
        col1, col2 = st.columns(2)
        with col1:
            # Get distinct source types
            source_types = [row[0] for row in db.query(Literature.source_type).distinct().all() if row[0]]
            source_type_filter = st.selectbox(
                "Filter by Source Type",
                ["All"] + source_types,
            )
        with col2:
            embedding_status_filter = st.selectbox(
                "Filter by Embedding Status",
                ["All", "Embedded", "Not Embedded", "In Progress"],
            )

        # Search
        search_term = st.text_input("Search by title, author, or DOI", "")

        # Query literature with eager loading of chunks
        query = db.query(Literature).options(joinedload(Literature.chunks))

        if source_type_filter != "All":
            query = query.filter(Literature.source_type == source_type_filter)

        if embedding_status_filter != "All":
            query = query.filter(Literature.embedding_status == embedding_status_filter)

        if search_term:
            query = query.filter(
                (Literature.title.ilike(f"%{search_term}%"))
                | (Literature.doi.ilike(f"%{search_term}%"))
                | (Literature.authors.any(lambda x: search_term.lower() in x.lower() if x else False))
            )

        literature_items = query.order_by(Literature.created_at.desc()).limit(100).all()

        st.metric("Literature Items Found", len(literature_items))
        if len(literature_items) == 100:
            st.info("Showing first 100 items. Use filters to narrow down.")

        # Access relationships while session is open
        literature_data_list = []
        for lit in literature_items:
            chunk_count = len(lit.chunks) if lit.chunks else 0
            literature_data_list.append(
                {
                    "lit": lit,
                    "chunk_count": chunk_count,
                }
            )

    # Display data outside session (relationships already loaded)
    if literature_data_list:
        # Summary table
        literature_data = []
        for item in literature_data_list:
            lit = item["lit"]
            literature_data.append(
                {
                    "Title": lit.title[:80] + "..." if len(lit.title) > 80 else lit.title,
                    "Source Type": lit.source_type or "-",
                    "Authors": ", ".join(lit.authors[:3]) if lit.authors else "-",
                    "Year": lit.year or "-",
                    "DOI": lit.doi or "-",
                    "Embedding Status": lit.embedding_status or "Not Embedded",
                    "Chunks": item["chunk_count"],
                    "Created": lit.created_at.strftime("%Y-%m-%d"),
                }
            )
        df_literature = pd.DataFrame(literature_data)
        st.dataframe(df_literature, width='stretch', hide_index=True)

        # Export button
        csv_literature = df_literature.to_csv(index=False)
        st.download_button(
            label="📥 Download Literature (CSV)",
            data=csv_literature,
            file_name="literature.csv",
            mime="text/csv",
        )

        st.markdown("---")
        st.subheader("Literature Details")

        for item in literature_data_list:
            lit = item["lit"]
            chunk_count = item["chunk_count"]

            with st.expander(f"**{lit.title}**"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{lit.id}`")
                    st.write(f"**Source Type:** {lit.source_type or '-'}")
                    if lit.authors:
                        st.write(f"**Authors:** {', '.join(lit.authors)}")
                    if lit.journal:
                        st.write(f"**Journal:** {lit.journal}")
                    if lit.year:
                        st.write(f"**Year:** {lit.year}")
                    if lit.doi:
                        st.write(f"**DOI:** {lit.doi}")
                with col2:
                    st.write(f"**Created:** {lit.created_at.strftime('%Y-%m-%d %H:%M')}")
                    st.write(f"**Embedding Status:** {lit.embedding_status or 'Not Embedded'}")
                    if lit.last_ingested_at:
                        st.write(f"**Last Ingested:** {lit.last_ingested_at.strftime('%Y-%m-%d %H:%M')}")
                    if lit.zotero_item_key:
                        st.write(f"**Zotero Key:** `{lit.zotero_item_key}`")

                if lit.abstract:
                    st.write("**Abstract:**")
                    st.text(lit.abstract[:500] + "..." if len(lit.abstract) > 500 else lit.abstract)

                if lit.tags:
                    st.write(f"**Tags:** {', '.join(lit.tags)}")

                # Linked chunks (already loaded via joinedload)
                if chunk_count > 0 and lit.chunks:
                    st.write(f"**RAG Chunks:** {chunk_count}")
                    # Show sample chunks
                    sample_chunks = lit.chunks[:5]
                    for chunk in sample_chunks:
                        with st.expander(
                            f"Chunk {chunk.chunk_index}: {chunk.snippet[:50]}..."
                            if chunk.snippet
                            else f"Chunk {chunk.chunk_index}"
                        ):
                            st.text(chunk.chunk_text[:500] + "..." if len(chunk.chunk_text) > 500 else chunk.chunk_text)
    else:
        st.info("No literature items found matching your filters.")
