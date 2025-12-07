"""Emails page for the Streamlit dashboard."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from sqlalchemy.orm import joinedload

from amprenta_rag.database.models import Email
from scripts.dashboard.db_session import db_session


def render_emails_page() -> None:
    """
    Render the Emails page with filtering, search, and detailed views.

    Features:
    - Filter by item type (Email, Note)
    - Search by title, sender
    - Display embedding status
    - Show linked RAG chunks
    - Display semantic metadata
    - Export functionality (CSV)
    """
    st.header("📧 Emails")

    with db_session() as db:
        # Filters
        col1, col2 = st.columns(2)
        with col1:
            # Get distinct item types
            item_types = [row[0] for row in db.query(Email.item_type).distinct().all() if row[0]]
            item_type_filter = st.selectbox(
                "Filter by Item Type",
                ["All"] + item_types,
            )
        with col2:
            embedding_status_filter = st.selectbox(
                "Filter by Embedding Status",
                ["All", "Embedded", "Not Embedded", "In Progress"],
            )

        # Search
        search_term = st.text_input("Search by title or sender", "")

        # Query emails with eager loading of chunks
        query = db.query(Email).options(joinedload(Email.chunks))

        if item_type_filter != "All":
            query = query.filter(Email.item_type == item_type_filter)

        if embedding_status_filter != "All":
            query = query.filter(Email.embedding_status == embedding_status_filter)

        if search_term:
            query = query.filter(
                (Email.title.ilike(f"%{search_term}%")) | (Email.from_sender.ilike(f"%{search_term}%"))
            )

        emails = query.order_by(Email.created_at.desc()).limit(100).all()

        st.metric("Emails Found", len(emails))
        if len(emails) == 100:
            st.info("Showing first 100 emails. Use filters to narrow down.")

        # Access relationships while session is open
        emails_data_list = []
        for email in emails:
            chunk_count = len(email.chunks) if email.chunks else 0
            emails_data_list.append(
                {
                    "email": email,
                    "chunk_count": chunk_count,
                }
            )

    # Display data outside session (relationships already loaded)
    if emails_data_list:
        # Summary table
        email_data = []
        for item in emails_data_list:
            email = item["email"]
            email_data.append(
                {
                    "Title": email.title[:60] + "..." if len(email.title) > 60 else email.title,
                    "From": email.from_sender or "-",
                    "Type": email.item_type or "-",
                    "Embedding Status": email.embedding_status or "Not Embedded",
                    "Chunks": item["chunk_count"],
                    "Created": email.created_at.strftime("%Y-%m-%d"),
                }
            )
        df_emails = pd.DataFrame(email_data)
        st.dataframe(df_emails, use_container_width=True, hide_index=True)

        # Export button
        csv_emails = df_emails.to_csv(index=False)
        st.download_button(
            label="📥 Download Emails (CSV)",
            data=csv_emails,
            file_name="emails.csv",
            mime="text/csv",
        )

        st.markdown("---")
        st.subheader("Email Details")

        for item in emails_data_list:
            email = item["email"]
            chunk_count = item["chunk_count"]

            with st.expander(f"**{email.title}**"):
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**ID:** `{email.id}`")
                    st.write(f"**Type:** {email.item_type or '-'}")
                    if email.from_sender:
                        st.write(f"**From:** {email.from_sender}")
                    if email.tags:
                        st.write(f"**Tags:** {', '.join(email.tags)}")
                with col2:
                    st.write(f"**Created:** {email.created_at.strftime('%Y-%m-%d %H:%M')}")
                    st.write(f"**Embedding Status:** {email.embedding_status or 'Not Embedded'}")
                    if email.last_ingested_at:
                        st.write(f"**Last Ingested:** {email.last_ingested_at.strftime('%Y-%m-%d %H:%M')}")

                # Semantic metadata
                if email.semantic_metadata:
                    st.write("**Semantic Metadata:**")
                    metadata = email.semantic_metadata
                    if isinstance(metadata, dict):
                        if metadata.get("diseases"):
                            st.caption(f"Diseases: {', '.join(metadata['diseases'])}")
                        if metadata.get("targets"):
                            st.caption(f"Targets: {', '.join(metadata['targets'])}")
                        if metadata.get("features"):
                            features = metadata["features"]
                            if isinstance(features, list):
                                st.caption(
                                    f"Features: {', '.join(features[:10])}"
                                    + (f" ... and {len(features) - 10} more" if len(features) > 10 else "")
                                )

                # Content preview
                if email.content:
                    st.write("**Content Preview:**")
                    st.text(email.content[:500] + "..." if len(email.content) > 500 else email.content)

                # Linked chunks (already loaded via joinedload)
                if chunk_count > 0 and email.chunks:
                    st.write(f"**RAG Chunks:** {chunk_count}")
                    # Show sample chunks
                    sample_chunks = email.chunks[:5]
                    for chunk in sample_chunks:
                        with st.expander(
                            f"Chunk {chunk.chunk_index}: {chunk.snippet[:50]}..."
                            if chunk.snippet
                            else f"Chunk {chunk.chunk_index}"
                        ):
                            st.text(chunk.chunk_text[:500] + "..." if len(chunk.chunk_text) > 500 else chunk.chunk_text)
    else:
        st.info("No emails found matching your filters.")
