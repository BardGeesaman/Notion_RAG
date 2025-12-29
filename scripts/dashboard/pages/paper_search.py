"""Scientific Paper Search and Ingestion page."""

from __future__ import annotations

import os
from typing import Any, Dict, List
from uuid import UUID

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 60) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_paper_search_page() -> None:
    """Render the Paper Search and Ingestion page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📄 Scientific Paper Search")
    st.caption("Search and ingest papers from PubMed, bioRxiv, and medRxiv.")

    tab1, tab2 = st.tabs(["Search Papers", "Ingested Papers"])

    with tab1:
        render_search_tab()

    with tab2:
        render_ingested_papers_tab()


def render_search_tab() -> None:
    """Render the paper search tab."""
    st.subheader("Search for Papers")

    col1, col2 = st.columns([3, 1])
    with col1:
        query = st.text_input(
            "Search Query",
            placeholder="Enter keywords, authors, or topics...",
            help="Search for scientific papers across PubMed and bioRxiv/medRxiv",
            key="paper_search_query",
        )
    with col2:
        source = st.selectbox(
            "Source",
            options=["pubmed", "biorxiv"],
            format_func=lambda x: "PubMed" if x == "pubmed" else "bioRxiv/medRxiv",
            help="Select the paper repository to search",
            key="paper_search_source",
        )

    col3, col4, col5 = st.columns([1, 1, 3])
    with col3:
        limit = st.number_input(
            "Max Results",
            min_value=1,
            max_value=100,
            value=10,
            step=1,
            help="Maximum number of search results to return",
            key="paper_search_limit",
        )
    with col4:
        search_btn = st.button(
            "🔍 Search",
            type="primary",
            use_container_width=True,
        )

    if search_btn and query:
        with st.spinner(f"Searching {source.upper()}..."):
            try:
                payload = {
                    "query": query,
                    "source": source,
                    "limit": limit,
                    "offset": 0,
                }
                response = _api_post("/api/v1/papers/search", payload)
                st.session_state["search_results"] = response.get("results", [])
                st.session_state["search_total"] = response.get("total", 0)
                st.success(f"Found {response.get('total', 0)} papers")
            except Exception as e:
                st.error(f"Search failed: {e}")
                st.session_state["search_results"] = []

    # Display search results
    if st.session_state.get("search_results"):
        results = st.session_state["search_results"]
        st.markdown(f"**{len(results)} Results**")

        for i, paper in enumerate(results):
            with st.container():
                st.markdown("---")
                col_title, col_action = st.columns([4, 1])

                with col_title:
                    title = paper.get("title", "No Title")
                    st.markdown(f"### {title}")

                    # Authors
                    authors = paper.get("authors", [])
                    if authors:
                        authors_str = ", ".join(authors[:3])
                        if len(authors) > 3:
                            authors_str += f" et al. ({len(authors)} authors)"
                        st.caption(f"**Authors:** {authors_str}")

                    # Journal and Year
                    journal = paper.get("journal", "")
                    year = paper.get("year", "")
                    if journal or year:
                        journal_year = f"**Journal:** {journal}" if journal else ""
                        if year:
                            journal_year += f" ({year})" if journal_year else f"**Year:** {year}"
                        st.caption(journal_year)

                    # Identifiers
                    pmid = paper.get("pmid", "")
                    doi = paper.get("doi", "")
                    id_str = []
                    if pmid:
                        id_str.append(f"PMID: {pmid}")
                    if doi:
                        id_str.append(f"DOI: {doi}")
                    if id_str:
                        st.caption(" | ".join(id_str))

                    # Abstract snippet
                    abstract = paper.get("abstract", "")
                    if abstract:
                        snippet = abstract[:300] + "..." if len(abstract) > 300 else abstract
                        with st.expander("Abstract", expanded=False):
                            st.write(abstract)

                with col_action:
                    st.write("")  # Spacer
                    st.write("")  # Spacer
                    ingest_btn = st.button(
                        "📥 Ingest",
                        key=f"ingest_btn_{i}",
                        use_container_width=True,
                        help="Ingest this paper into the database",
                    )

                    if ingest_btn:
                        with st.spinner("Ingesting..."):
                            try:
                                ingest_payload = {
                                    "pmid": paper.get("pmid"),
                                    "doi": paper.get("doi"),
                                    "fetch_fulltext": False,  # TODO: Enable after PMC integration
                                }
                                result = _api_post("/api/v1/papers/ingest", ingest_payload)
                                if result.get("already_exists"):
                                    st.info("Paper already ingested")
                                else:
                                    st.success(f"Ingested! ID: {result.get('literature_id')}")
                                    # Store ingested ID for viewing
                                    if "ingested_ids" not in st.session_state:
                                        st.session_state["ingested_ids"] = []
                                    st.session_state["ingested_ids"].append(result.get("literature_id"))
                            except Exception as e:
                                st.error(f"Ingestion failed: {e}")


def render_ingested_papers_tab() -> None:
    """Render the ingested papers tab."""
    st.subheader("View Ingested Papers")

    # Input for literature ID
    lit_id_input = st.text_input(
        "Literature ID",
        placeholder="Enter UUID of ingested paper...",
        help="UUID of a paper previously ingested",
        key="ingested_paper_id",
    )

    view_btn = st.button("📖 View Paper", type="primary")

    if view_btn and lit_id_input:
        try:
            literature_id = UUID(lit_id_input)
        except ValueError:
            st.error("Invalid UUID format")
            return

        with st.spinner("Fetching paper details..."):
            try:
                paper = _api_get(f"/api/v1/papers/{literature_id}")
                st.session_state["viewed_paper"] = paper
            except Exception as e:
                st.error(f"Failed to fetch paper: {e}")
                return

    # Display viewed paper
    if st.session_state.get("viewed_paper"):
        paper = st.session_state["viewed_paper"]

        st.markdown(f"## {paper.get('title', 'No Title')}")

        # Metadata
        col1, col2, col3 = st.columns(3)
        with col1:
            if paper.get("pmid"):
                st.metric("PMID", paper["pmid"])
        with col2:
            if paper.get("source"):
                st.metric("Source", paper["source"].upper())
        with col3:
            if paper.get("full_text_available"):
                st.metric("Full Text", "✅ Available" if paper["full_text_available"] else "❌ Not Available")

        # DOI and URL
        if paper.get("doi"):
            st.caption(f"**DOI:** {paper['doi']}")
        if paper.get("url"):
            st.caption(f"**URL:** [{paper['url']}]({paper['url']})")

        # MeSH Terms
        mesh_terms = paper.get("mesh_terms", [])
        if mesh_terms:
            st.caption(f"**MeSH Terms:** {', '.join(mesh_terms)}")

        # Abstract
        if paper.get("abstract"):
            st.markdown("### Abstract")
            st.write(paper["abstract"])

        # Sections (if full text available)
        sections = paper.get("sections", [])
        if sections:
            st.markdown("### Full Text Sections")
            for section in sections:
                with st.expander(f"📄 {section['title']}", expanded=False):
                    st.write(section["content"])
        else:
            if paper.get("full_text_available"):
                st.info("Full text is available but sections are not yet parsed.")
            else:
                st.info("Full text is not available for this paper. Only abstract is shown.")
        
        # Extraction Actions
        st.markdown("---")
        st.markdown("### Data Extraction")
        
        col1, col2 = st.columns(2)
        with col1:
            if st.button("🔬 Extract Experiments", help="Extract structured experiment data from PDF"):
                with st.spinner("Extracting experiments..."):
                    try:
                        result = _api_post(f"/api/v1/papers/{literature_id}/extract", {})
                        st.success(f"Extracted {result.get('experiments_extracted', 0)} experiments")
                        # Refresh experiments list
                        st.session_state.pop("viewed_experiments", None)
                    except Exception as e:
                        st.error(f"Extraction failed: {e}")
        
        with col2:
            if st.button("📊 View Extracted Experiments"):
                with st.spinner("Loading extracted experiments..."):
                    try:
                        experiments = _api_get(f"/api/v1/papers/{literature_id}/experiments")
                        st.session_state["viewed_experiments"] = experiments
                    except Exception as e:
                        st.error(f"Failed to load experiments: {e}")
        
        # Display extracted experiments
        if st.session_state.get("viewed_experiments"):
            experiments = st.session_state["viewed_experiments"]
            
            if experiments:
                st.markdown("#### Extracted Experiments")
                for i, exp in enumerate(experiments, 1):
                    with st.expander(f"Experiment {i}: {exp.get('experiment_type', 'Unknown')}"):
                        if exp.get("cell_line"):
                            st.markdown(f"**Cell Line:** {exp['cell_line']}")
                        if exp.get("treatment"):
                            st.markdown(f"**Treatment:** {exp['treatment']}")
                        if exp.get("concentration"):
                            st.markdown(f"**Concentration:** {exp['concentration']}")
                        if exp.get("timepoint"):
                            st.markdown(f"**Timepoint:** {exp['timepoint']}")
                        if exp.get("measured_entities"):
                            st.markdown(f"**Measured:** {', '.join(exp['measured_entities'][:10])}")
                        if exp.get("key_findings"):
                            st.markdown(f"**Findings:** {exp['key_findings']}")
                        if exp.get("extraction_confidence"):
                            st.progress(exp["extraction_confidence"] / 100.0, text=f"Confidence: {exp['extraction_confidence']}%")
            else:
                st.info("No experiments extracted yet. Click 'Extract Experiments' to analyze this paper.")

    # Show recently ingested papers
    if st.session_state.get("ingested_ids"):
        st.markdown("---")
        st.subheader("Recently Ingested Papers")
        for lit_id in st.session_state["ingested_ids"][-5:]:
            st.code(lit_id, language=None)


if __name__ == "__main__":
    render_paper_search_page()

