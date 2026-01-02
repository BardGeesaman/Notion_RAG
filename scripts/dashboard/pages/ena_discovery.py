"""ENA Discovery Dashboard for browsing and ingesting genomics studies."""

from __future__ import annotations

import streamlit as st

from scripts.dashboard.core.api import _api_get, _api_post
from scripts.dashboard.core.state import require_auth


def render_ena_discovery_page():
    """Main entry point for ENA Discovery page."""
    require_auth()
    
    st.title("🧬 ENA Discovery")
    st.caption("Search and browse European Nucleotide Archive for genomics studies")
    
    # Tab layout
    tab_search, tab_results, tab_ingest = st.tabs([
        "🔍 Search",
        "📋 Results",
        "⬇️ Ingest"
    ])
    
    with tab_search:
        _render_search_tab()
    
    with tab_results:
        _render_results_tab()
    
    with tab_ingest:
        _render_ingest_tab()


def _render_search_tab():
    """Search tab for ENA queries."""
    st.subheader("🔍 Search ENA")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        query = st.text_input(
            "Search Keywords",
            placeholder="e.g., Homo sapiens, cancer, RNA-Seq",
            key="ena_query"
        )
    
    with col2:
        max_results = st.slider("Max Results", 5, 100, 20, key="ena_max_results")
    
    col3, col4 = st.columns(2)
    
    with col3:
        organism = st.text_input(
            "Organism Filter",
            placeholder="e.g., Homo sapiens",
            key="ena_organism"
        )
    
    with col4:
        library_strategy = st.selectbox(
            "Library Strategy",
            ["", "RNA-Seq", "WGS", "WXS", "ChIP-Seq", "ATAC-seq", "Hi-C"],
            key="ena_library"
        )
    
    if query and st.button("🔍 Search ENA", key="ena_search_btn", type="primary"):
        with st.spinner("Searching ENA..."):
            params = f"?q={query}&max_results={max_results}"
            if organism:
                params += f"&organism={organism}"
            if library_strategy:
                params += f"&library_strategy={library_strategy}"
            
            response = _api_get(f"/genomics/ena/search{params}")
            
            if response and response.get("success"):
                st.session_state["ena_search_results"] = response.get("results", [])
                st.success(f"✅ Found {response.get('total', 0)} studies")
            else:
                st.error(f"❌ Search failed: {response.get('error', 'Unknown error') if response else 'API error'}")


def _render_results_tab():
    """Results tab showing search results."""
    st.subheader("📋 Search Results")
    
    results = st.session_state.get("ena_search_results", [])
    
    if not results:
        st.info("No results yet. Use the Search tab to find studies.")
        return
    
    st.metric("Studies Found", len(results))
    
    # Initialize selected studies in session state
    if "ena_selected_studies" not in st.session_state:
        st.session_state["ena_selected_studies"] = []
    
    for i, study in enumerate(results):
        with st.expander(f"📁 {study.get('run_accession', 'Unknown')} - {study.get('title', 'No title')[:60]}..."):
            col1, col2 = st.columns([3, 1])
            
            with col1:
                st.write(f"**Run Accession:** {study.get('run_accession', 'N/A')}")
                st.write(f"**Study Accession:** {study.get('study_accession', 'N/A')}")
                st.write(f"**Organism:** {study.get('organism', 'N/A')}")
                st.write(f"**Platform:** {study.get('platform', 'N/A')}")
                st.write(f"**FASTQ Files:** {study.get('fastq_count', 0)}")
            
            with col2:
                selected = st.checkbox(
                    "Select for Ingest",
                    key=f"ena_select_{i}",
                    value=study.get('run_accession') in st.session_state["ena_selected_studies"]
                )
                
                if selected and study.get('run_accession') not in st.session_state["ena_selected_studies"]:
                    st.session_state["ena_selected_studies"].append(study.get('run_accession'))
                elif not selected and study.get('run_accession') in st.session_state["ena_selected_studies"]:
                    st.session_state["ena_selected_studies"].remove(study.get('run_accession'))


def _render_ingest_tab():
    """Ingest tab for triggering study ingestion."""
    st.subheader("⬇️ Ingest Selected Studies")
    
    selected = st.session_state.get("ena_selected_studies", [])
    
    if not selected:
        st.info("No studies selected. Use the Results tab to select studies for ingestion.")
        return
    
    st.metric("Selected Studies", len(selected))
    
    st.write("**Studies to ingest:**")
    for study_id in selected:
        st.write(f"- {study_id}")
    
    st.warning("⚠️ FASTQ files can be very large (GB to TB). Ingestion may take significant time.")
    
    if st.button("🚀 Start Ingestion", key="ena_ingest_btn", type="primary"):
        with st.spinner("Queuing studies for ingestion..."):
            response = _api_post("/genomics/ena/ingest", {"study_ids": selected})
            
            if response and response.get("success"):
                st.success(f"✅ {response.get('message', 'Ingestion queued')}")
                # Clear selection after successful queue
                st.session_state["ena_selected_studies"] = []
            else:
                st.error(f"❌ Ingestion failed: {response.get('error', 'Unknown error') if response else 'API error'}")
