"""Ontology Management page for managing controlled vocabularies."""

from __future__ import annotations

from typing import List

import pandas as pd
import streamlit as st

from amprenta_rag.database.models import OntologyTerm
from amprenta_rag.auth.session import get_current_user
from scripts.dashboard.db_session import db_session


def render_ontology_page() -> None:
    """Render the Ontology Management page."""
    user = get_current_user()

    # Admin only
    if not user or user.get("role") != "admin":
        st.error("Access denied. Only administrators can manage ontologies.")
        return

    st.header("📚 Ontology Management")
    st.markdown("Manage controlled vocabularies and ontology terms.")

    tabs = st.tabs(["Browse", "Add Term"])

    with tabs[0]:
        _render_browse_tab()

    with tabs[1]:
        _render_add_term_tab()


def _render_browse_tab() -> None:
    """Render the browse tab with vocabulary filter and term display."""
    with db_session() as db:
        # Get unique vocabularies
        vocabularies = db.query(OntologyTerm.vocabulary).distinct().all()
        vocabulary_options = ["All"] + [v[0] for v in vocabularies if v[0]]

        selected_vocabulary = st.selectbox("Filter by Vocabulary", vocabulary_options)

        # Build query
        query = db.query(OntologyTerm)

        if selected_vocabulary != "All":
            query = query.filter(OntologyTerm.vocabulary == selected_vocabulary)

        # Filter active only by default
        show_inactive = st.checkbox("Show Inactive Terms", value=False)
        if not show_inactive:
            query = query.filter(OntologyTerm.is_active)

        terms = query.order_by(OntologyTerm.vocabulary, OntologyTerm.term).all()

        st.metric("Total Terms", len(terms))

        if terms:
            # Display as table
            term_data = []
            for term in terms:
                parent_term = ""
                if term.parent_id:
                    parent = db.query(OntologyTerm).filter(OntologyTerm.id == term.parent_id).first()
                    parent_term = parent.term if parent else "Unknown"

                term_data.append({
                    "Vocabulary": term.vocabulary,
                    "Term": term.term,
                    "Description": term.description[:100] + "..." if term.description and len(term.description) > 100 else (term.description or ""),
                    "Parent": parent_term or "None",
                    "Status": "✅ Active" if term.is_active else "❌ Inactive",
                    "Created": term.created_at.strftime("%Y-%m-%d") if term.created_at else "",
                })

            df_terms = pd.DataFrame(term_data)
            st.dataframe(df_terms, use_container_width=True, hide_index=True)

            # Tree view option
            st.markdown("---")
            if st.checkbox("Show Tree View", value=False):
                st.subheader("Hierarchical Tree View")
                _render_term_tree(terms, db)
        else:
            st.info("No terms found. Add terms using the 'Add Term' tab.")


def _render_term_tree(terms: List[OntologyTerm], db) -> None:
    """Render terms in a hierarchical tree structure."""
    # Group by vocabulary
    vocab_groups = {}
    for term in terms:
        if term.vocabulary not in vocab_groups:
            vocab_groups[term.vocabulary] = []
        vocab_groups[term.vocabulary].append(term)

    for vocabulary, vocab_terms in vocab_groups.items():
        st.markdown(f"### {vocabulary}")

        # Find root terms (no parent)
        root_terms = [t for t in vocab_terms if not t.parent_id]

        # Build tree structure
        def render_term_recursive(term: OntologyTerm, level: int = 0):
            indent = "  " * level
            status_icon = "✅" if term.is_active else "❌"
            st.markdown(f"{indent}- {status_icon} **{term.term}**")
            if term.description:
                st.caption(f"{indent}  {term.description[:100]}{'...' if len(term.description) > 100 else ''}")

            # Find children
            children = [t for t in vocab_terms if t.parent_id == term.id]
            for child in children:
                render_term_recursive(child, level + 1)

        for root_term in root_terms:
            render_term_recursive(root_term)


def _render_add_term_tab() -> None:
    """Render the add term form."""
    with db_session() as db:
        st.subheader("Add New Ontology Term")

        with st.form("add_ontology_term", clear_on_submit=True):
            vocabulary = st.text_input("Vocabulary*", placeholder="e.g., disease, tissue, cell_type")

            term = st.text_input("Term*", placeholder="e.g., Alzheimer's Disease, Brain, Neuron")

            description = st.text_area("Description (optional)", placeholder="Description of this term...", height=100)

            # Parent selector
            link_parent = st.checkbox("Link to Parent Term", value=False)
            parent_id = None

            if link_parent:
                # Get existing terms for parent selection
                existing_terms = db.query(OntologyTerm).filter(OntologyTerm.is_active).order_by(OntologyTerm.vocabulary, OntologyTerm.term).all()
                if existing_terms:
                    parent_options = {f"{t.vocabulary}: {t.term}": t.id for t in existing_terms}
                    selected_parent_label = st.selectbox("Select Parent Term", list(parent_options.keys()))
                    parent_id = parent_options[selected_parent_label]
                else:
                    st.info("No existing terms available for parent selection.")

            is_active = st.checkbox("Active", value=True)

            submitted = st.form_submit_button("💾 Save Term", type="primary")

            if submitted:
                if not vocabulary or not vocabulary.strip():
                    st.error("Vocabulary is required.")
                elif not term or not term.strip():
                    st.error("Term is required.")
                else:
                    try:
                        # Check for duplicate
                        existing = db.query(OntologyTerm).filter(
                            OntologyTerm.vocabulary == vocabulary.strip(),
                            OntologyTerm.term == term.strip()
                        ).first()

                        if existing:
                            st.error(f"Term '{term}' already exists in vocabulary '{vocabulary}'.")
                        else:
                            ontology_term = OntologyTerm(
                                vocabulary=vocabulary.strip(),
                                term=term.strip(),
                                description=description.strip() if description.strip() else None,
                                parent_id=parent_id,
                                is_active=is_active,
                            )

                            db.add(ontology_term)
                            db.commit()
                            st.success(f"Term '{term}' added successfully!")
                            st.rerun()

                    except Exception as e:
                        st.error(f"Failed to add term: {e}")


__all__ = ["render_ontology_page"]

