"""Compound Scoring (Relevance & Novelty) page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 120) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_scoring_page() -> None:
    """Render the Compound Scoring page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("📊 Compound Scoring")
    st.caption("Score compounds for relevance and novelty using AI-powered analysis.")
    
    tab1, tab2, tab3 = st.tabs(["Single Compound", "Batch Scoring", "Results"])
    
    with tab1:
        render_single_scoring_tab()
    
    with tab2:
        render_batch_scoring_tab()
    
    with tab3:
        render_results_tab()


def render_single_scoring_tab() -> None:
    """Render single compound scoring tab."""
    st.subheader("Score Single Compound")
    
    # Compound details
    st.markdown("### Compound Details")
    compound_id = st.text_input("Compound ID", placeholder="e.g., CMPD001")
    
    col1, col2 = st.columns(2)
    with col1:
        title = st.text_input("Title/Name", placeholder="e.g., Novel ALS Compound")
        species = st.text_input("Species", placeholder="e.g., human")
    with col2:
        assay_type = st.text_input("Assay Type", placeholder="e.g., cell_viability")
        sample_count = st.number_input("Sample Count", min_value=0, value=50, step=1)
    
    description = st.text_area(
        "Description",
        height=100,
        placeholder="Describe the compound or dataset...",
    )
    
    # Research context
    st.markdown("### Research Context (for Relevance)")
    
    diseases = st.text_input("Diseases (comma-separated)", placeholder="e.g., ALS, Parkinson's")
    targets = st.text_input("Targets (comma-separated)", placeholder="e.g., SOD1, TDP43")
    
    col1, col2 = st.columns(2)
    with col1:
        context_species = st.text_input("Context Species", placeholder="e.g., human, mouse", key="context_species")
    with col2:
        min_sample_size = st.number_input("Min Sample Size", min_value=1, value=20, step=5)
    
    # Scoring options
    st.markdown("### Scoring Options")
    col1, col2 = st.columns(2)
    with col1:
        score_relevance = st.checkbox("Score Relevance", value=True)
    with col2:
        score_novelty = st.checkbox("Score Novelty", value=True)
    
    if st.button("Score Compound", type="primary", disabled=not compound_id):
        results = {}
        
        # Relevance scoring
        if score_relevance:
            with st.spinner("Scoring relevance..."):
                try:
                    relevance_result = _api_post(
                        "/api/v1/score/relevance",
                        {
                            "item": {
                                "id": compound_id,
                                "title": title or compound_id,
                                "description": description,
                                "species": species,
                                "assay_type": assay_type,
                                "sample_count": sample_count,
                            },
                            "context": {
                                "diseases": [d.strip() for d in diseases.split(",")] if diseases else [],
                                "targets": [t.strip() for t in targets.split(",")] if targets else [],
                                "species": [s.strip() for s in context_species.split(",")] if context_species else [],
                                "assay_types": [],
                                "min_sample_size": min_sample_size,
                            },
                            "criteria": None,
                        },
                        timeout=60,
                    )
                    results["relevance"] = relevance_result
                except Exception as e:
                    st.error(f"Relevance scoring failed: {e}")
        
        # Novelty scoring
        if score_novelty:
            with st.spinner("Scoring novelty..."):
                try:
                    novelty_result = _api_post(
                        "/api/v1/score/novelty",
                        {
                            "item": {
                                "id": compound_id,
                                "title": title or compound_id,
                                "description": description,
                                "species": species,
                                "assay_type": assay_type,
                                "sample_count": sample_count,
                            },
                            "existing_items": [],  # Would load from DB in production
                        },
                        timeout=60,
                    )
                    results["novelty"] = novelty_result
                except Exception as e:
                    st.error(f"Novelty scoring failed: {e}")
        
        if results:
            st.session_state["scoring_result"] = results
    
    # Display results
    result = st.session_state.get("scoring_result")
    if result:
        st.markdown("---")
        st.subheader("Scoring Results")
        
        if "relevance" in result:
            st.markdown("### Relevance Score")
            rel = result["relevance"]
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Overall", f"{rel.get('overall_score', 0.0):.2f}")
            with col2:
                st.metric("Disease Match", f"{rel.get('disease_match', 0.0):.2f}")
            with col3:
                st.metric("Target Overlap", f"{rel.get('target_overlap', 0.0):.2f}")
            with col4:
                st.metric("Data Quality", f"{rel.get('data_quality', 0.0):.2f}")
            
            if rel.get("explanation"):
                with st.expander("Explanation"):
                    st.write(rel["explanation"])
        
        if "novelty" in result:
            st.markdown("### Novelty Score")
            nov = result["novelty"]
            
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Novelty Score", f"{nov.get('novelty_score', 0.0):.2f}")
            with col2:
                st.metric("Max Similarity", f"{nov.get('max_similarity', 0.0):.2f}")
            
            if nov.get("explanation"):
                with st.expander("Explanation"):
                    st.write(nov["explanation"])


def render_batch_scoring_tab() -> None:
    """Render batch scoring tab."""
    st.subheader("Batch Scoring")
    
    st.markdown("Score multiple compounds at once for ranking and prioritization.")
    
    # SMILES input
    smiles_input = st.text_area(
        "Compound SMILES (one per line)",
        height=200,
        placeholder="CCO\nCC(C)O\nC1CCCCC1",
        help="Enter SMILES strings, one per line",
    )
    
    # Context
    col1, col2 = st.columns(2)
    with col1:
        diseases = st.text_input("Diseases", placeholder="ALS, Parkinson's", key="batch_diseases")
        targets = st.text_input("Targets", placeholder="SOD1, TDP43", key="batch_targets")
    with col2:
        species = st.text_input("Species", placeholder="human, mouse", key="batch_species")
        min_samples = st.number_input("Min Samples", min_value=1, value=20, key="batch_min")
    
    # Options
    col1, col2 = st.columns(2)
    with col1:
        score_relevance = st.checkbox("Score Relevance", value=True, key="batch_rel")
    with col2:
        score_novelty = st.checkbox("Score Novelty", value=True, key="batch_nov")
    
    if st.button("Score Batch", type="primary", disabled=not smiles_input):
        smiles_list = [s.strip() for s in smiles_input.strip().split("\n") if s.strip()]
        
        with st.spinner(f"Scoring {len(smiles_list)} compounds..."):
            try:
                # Build items list
                items = [
                    {
                        "id": f"compound_{i}",
                        "title": f"Compound {i}",
                        "description": smiles,
                        "species": species or None,
                        "assay_type": "unknown",
                        "sample_count": min_samples,
                    }
                    for i, smiles in enumerate(smiles_list, 1)
                ]
                
                result = _api_post(
                    "/api/v1/score/batch",
                    {
                        "items": items,
                        "context": {
                            "diseases": [d.strip() for d in diseases.split(",")] if diseases else [],
                            "targets": [t.strip() for t in targets.split(",")] if targets else [],
                            "species": [s.strip() for s in species.split(",")] if species else [],
                            "assay_types": [],
                            "min_sample_size": min_samples,
                        },
                        "score_relevance": score_relevance,
                        "score_novelty": score_novelty,
                    },
                    timeout=180,
                )
                st.session_state["batch_scoring_result"] = result
            except Exception as e:
                st.error(f"Batch scoring failed: {e}")
                return
    
    # Display results
    result = st.session_state.get("batch_scoring_result")
    if result:
        st.success(f"Scored {result.get('total_items', 0)} compounds in {result.get('processing_time_seconds', 0.0):.1f}s")
        
        items = result.get("items", [])
        if items:
            # Build results table
            rows = []
            for item in items:
                row = {"Compound ID": item.get("item_id", "")}
                
                if score_relevance and item.get("relevance_score"):
                    rel = item["relevance_score"]
                    row["Relevance"] = f"{rel.get('overall_score', 0.0):.2f}"
                    row["Disease Match"] = f"{rel.get('disease_match', 0.0):.2f}"
                
                if score_novelty and item.get("novelty_score"):
                    nov = item["novelty_score"]
                    row["Novelty"] = f"{nov.get('novelty_score', 0.0):.2f}"
                
                rows.append(row)
            
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # Download
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 Download Scores (CSV)",
                data=csv,
                file_name="compound_scores.csv",
                mime="text/csv",
            )


def render_results_tab() -> None:
    """Render results history tab."""
    st.subheader("Scoring Results History")
    
    # Show recent scoring results from session state
    single_result = st.session_state.get("scoring_result")
    batch_result = st.session_state.get("batch_scoring_result")
    
    if single_result:
        st.markdown("### Recent Single Compound Score")
        if "relevance" in single_result:
            st.metric("Relevance", f"{single_result['relevance'].get('overall_score', 0.0):.2f}")
        if "novelty" in single_result:
            st.metric("Novelty", f"{single_result['novelty'].get('novelty_score', 0.0):.2f}")
    
    if batch_result:
        st.markdown("### Recent Batch Scoring")
        st.markdown(f"**Total Items**: {batch_result.get('total_items', 0)}")
        st.markdown(f"**Processing Time**: {batch_result.get('processing_time_seconds', 0.0):.1f}s")
    
    if not single_result and not batch_result:
        st.info("No scoring results yet. Run scoring from the other tabs.")


if __name__ == "__main__":
    render_scoring_page()

