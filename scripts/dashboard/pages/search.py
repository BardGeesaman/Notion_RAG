"""Global Search page for the Streamlit dashboard."""

from __future__ import annotations

from typing import Any, Dict, List

import streamlit as st

from amprenta_rag.database.models import Dataset, Email, Experiment, Feature, Literature, Program, Signature
from scripts.dashboard.db_session import db_session


def render_search_page() -> None:
    """
    Render the Global Search page.

    Features:
    - Search across all entity types
    - Filter by entity type
    - Highlight search terms
    - Quick navigation to entity details
    """
    st.header("🔍 Global Search")
    st.markdown(
        "Search across all entities in the system: datasets, programs, experiments, features, signatures, literature, and emails."
    )

    # Search input
    search_query = st.text_input(
        "Search",
        placeholder="Enter search term...",
        key="global_search",
    )

    # Entity type filter
    entity_types = st.multiselect(
        "Filter by entity type",
        [
            "Datasets",
            "Programs",
            "Experiments",
            "Features",
            "Signatures",
            "Literature",
            "Emails",
        ],
        default=["Datasets", "Programs", "Experiments", "Features", "Signatures"],
    )

    if not search_query:
        st.info("Enter a search term to find entities across the system.")
        return

    if not entity_types:
        st.warning("Please select at least one entity type to search.")
        return

    # Perform search
    with st.spinner("Searching..."):
        results = perform_global_search(search_query, entity_types)

    # Display results
    if not results:
        st.info("No results found. Try a different search term or expand entity type filters.")
        return

    # Summary
    total_results = sum(len(entity_results) for entity_results in results.values())
    st.success(f"Found {total_results} result(s)")

    # Display results by entity type
    for entity_type, entity_results in results.items():
        if not entity_results:
            continue

        with st.expander(f"{entity_type} ({len(entity_results)})", expanded=True):
            for result in entity_results[:50]:  # Limit to 50 per type
                display_search_result(entity_type, result)

            if len(entity_results) > 50:
                st.caption(f"... and {len(entity_results) - 50} more {entity_type.lower()}")


def perform_global_search(
    query: str,
    entity_types: List[str],
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Perform global search across selected entity types.

    Args:
        query: Search query string
        entity_types: List of entity types to search

    Returns:
        Dictionary mapping entity type to list of results
    """
    results: Dict[str, List[Dict[str, Any]]] = {}

    with db_session() as db:
        # Search Datasets
        if "Datasets" in entity_types:
            datasets = (
                db.query(Dataset)
                .filter(
                    (Dataset.name.ilike(f"%{query}%"))
                    | (Dataset.description.ilike(f"%{query}%"))
                    | (Dataset.omics_type.ilike(f"%{query}%"))
                )
                .limit(100)
                .all()
            )
            results["Datasets"] = [
                {
                    "id": str(ds.id),
                    "name": ds.name,
                    "type": ds.omics_type,
                    "description": ds.description or "",
                    "created": ds.created_at.strftime("%Y-%m-%d") if ds.created_at else "",
                }
                for ds in datasets
            ]

        # Search Programs
        if "Programs" in entity_types:
            programs = (
                db.query(Program)
                .filter((Program.name.ilike(f"%{query}%")) | (Program.description.ilike(f"%{query}%")))
                .limit(100)
                .all()
            )
            results["Programs"] = [
                {
                    "id": str(p.id),
                    "name": p.name,
                    "description": p.description or "",
                    "created": p.created_at.strftime("%Y-%m-%d") if p.created_at else "",
                }
                for p in programs
            ]

        # Search Experiments
        if "Experiments" in entity_types:
            experiments = (
                db.query(Experiment)
                .filter(
                    (Experiment.name.ilike(f"%{query}%"))
                    | (Experiment.description.ilike(f"%{query}%"))
                    | (Experiment.type.ilike(f"%{query}%"))
                )
                .limit(100)
                .all()
            )
            results["Experiments"] = [
                {
                    "id": str(e.id),
                    "name": e.name,
                    "type": e.type or "",
                    "description": e.description or "",
                    "created": e.created_at.strftime("%Y-%m-%d") if e.created_at else "",
                }
                for e in experiments
            ]

        # Search Features
        if "Features" in entity_types:
            features = (
                db.query(Feature)
                .filter((Feature.name.ilike(f"%{query}%")) | (Feature.feature_type.ilike(f"%{query}%")))
                .limit(100)
                .all()
            )
            results["Features"] = [
                {
                    "id": str(f.id),
                    "name": f.name,
                    "type": f.feature_type,
                    "created": f.created_at.strftime("%Y-%m-%d") if f.created_at else "",
                }
                for f in features
            ]

        # Search Signatures
        if "Signatures" in entity_types:
            signatures = (
                db.query(Signature)
                .filter((Signature.name.ilike(f"%{query}%")) | (Signature.description.ilike(f"%{query}%")))
                .limit(100)
                .all()
            )
            results["Signatures"] = [
                {
                    "id": str(s.id),
                    "name": s.name,
                    "description": s.description or "",
                    "created": s.created_at.strftime("%Y-%m-%d") if s.created_at else "",
                }
                for s in signatures
            ]

        # Search Literature
        if "Literature" in entity_types:
            literature = (
                db.query(Literature)
                .filter(
                    (Literature.title.ilike(f"%{query}%"))
                    | (Literature.abstract.ilike(f"%{query}%"))
                    | (Literature.journal.ilike(f"%{query}%"))
                )
                .limit(100)
                .all()
            )
            results["Literature"] = [
                {
                    "id": str(l.id),
                    "name": l.title or "Untitled",
                    "journal": l.journal or "",
                    "year": l.year or "",
                    "created": l.created_at.strftime("%Y-%m-%d") if l.created_at else "",
                }
                for l in literature
            ]

        # Search Emails
        if "Emails" in entity_types:
            emails = (
                db.query(Email)
                .filter((Email.title.ilike(f"%{query}%")) | (Email.from_sender.ilike(f"%{query}%")))
                .limit(100)
                .all()
            )
            results["Emails"] = [
                {
                    "id": str(e.id),
                    "name": e.title or "Untitled",
                    "from": e.from_sender or "",
                    "created": e.created_at.strftime("%Y-%m-%d") if e.created_at else "",
                }
                for e in emails
            ]

    return results


def display_search_result(_entity_type: str, result: Dict[str, Any]) -> None:
    """Display a single search result."""
    col1, col2 = st.columns([3, 1])

    with col1:
        st.write(f"**{result['name']}**")
        if result.get("description"):
            st.caption(
                result["description"][:200] + "..."
                if len(result.get("description", "")) > 200
                else result["description"]
            )
        if result.get("type"):
            st.caption(f"Type: {result['type']}")
        if result.get("journal"):
            st.caption(f"Journal: {result['journal']}")
        if result.get("from"):
            st.caption(f"From: {result['from']}")

    with col2:
        st.caption(f"ID: `{result['id'][:8]}...`")
        if result.get("created"):
            st.caption(f"Created: {result['created']}")

    st.markdown("---")
