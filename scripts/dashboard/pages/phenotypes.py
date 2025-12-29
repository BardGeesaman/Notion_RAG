"""Phenotype (HPO) Exploration page."""

from __future__ import annotations

import os
from typing import Any, Dict, List

import httpx
import pandas as pd
import streamlit as st


API_BASE = os.environ.get("API_URL", "http://localhost:8000")


def _api_get(path: str, *, timeout: int = 30) -> Any:
    """Make GET request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.get(f"{API_BASE}{path}")
    r.raise_for_status()
    return r.json()


def _api_post(path: str, payload: Dict[str, Any], *, timeout: int = 30) -> Any:
    """Make POST request to API."""
    with httpx.Client(timeout=timeout) as client:
        r = client.post(f"{API_BASE}{path}", json=payload)
    r.raise_for_status()
    return r.json()


def render_phenotypes_page() -> None:
    """Render the Phenotypes Exploration page."""
    from scripts.dashboard.auth import require_auth
    require_auth()
    
    st.header("🧬 Phenotype Explorer")
    st.caption("Explore Human Phenotype Ontology (HPO) terms and gene-phenotype associations.")
    
    tab1, tab2 = st.tabs(["HPO Term Lookup", "Query Expansion"])
    
    with tab1:
        render_hpo_lookup_tab()
    
    with tab2:
        render_query_expansion_tab()


def render_hpo_lookup_tab() -> None:
    """Render HPO term lookup tab."""
    st.subheader("HPO Term → Genes")
    
    st.markdown("Enter an HPO term ID to find associated genes.")
    
    hpo_id = st.text_input(
        "HPO Term ID",
        placeholder="e.g., HP:0000001",
        help="Human Phenotype Ontology term identifier",
    )
    
    if st.button("Lookup Genes", type="primary", disabled=not hpo_id):
        with st.spinner(f"Looking up genes for {hpo_id}..."):
            try:
                genes = _api_get(f"/api/phenotypes/{hpo_id}/genes")
                st.session_state["hpo_lookup_result"] = {
                    "hpo_id": hpo_id,
                    "genes": genes,
                }
            except httpx.HTTPStatusError as e:
                if e.response.status_code == 404:
                    st.error(f"HPO term '{hpo_id}' not found or no gene associations.")
                else:
                    st.error(f"Lookup failed: {e}")
                return
            except Exception as e:
                st.error(f"Lookup failed: {e}")
                return
    
    # Display results
    result = st.session_state.get("hpo_lookup_result")
    if result:
        hpo_id = result.get("hpo_id", "")
        genes = result.get("genes", [])
        
        st.success(f"Found {len(genes)} genes associated with {hpo_id}")
        
        if genes:
            # Display as columns
            num_cols = 4
            cols = st.columns(num_cols)
            for i, gene in enumerate(genes):
                with cols[i % num_cols]:
                    st.markdown(f"- {gene}")
            
            # Download as CSV
            df = pd.DataFrame({"Gene": genes})
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 Download Genes (CSV)",
                data=csv,
                file_name=f"genes_{hpo_id}.csv",
                mime="text/csv",
            )
        else:
            st.info("No genes found for this HPO term.")


def render_query_expansion_tab() -> None:
    """Render query expansion tab."""
    st.subheader("Clinical Query → HPO Terms + Genes")
    
    st.markdown("Enter a clinical description to extract HPO terms and find associated genes.")
    
    query = st.text_area(
        "Clinical Query",
        height=150,
        placeholder="e.g., Patient has muscle weakness and progressive motor neuron degeneration...",
        help="Describe phenotypes in natural language",
    )
    
    if st.button("Expand Query", type="primary", disabled=not query):
        with st.spinner("Expanding query..."):
            try:
                result = _api_post("/api/phenotypes/expand-query", {"query": query})
                st.session_state["query_expansion_result"] = result
            except Exception as e:
                st.error(f"Query expansion failed: {e}")
                return
    
    # Display results
    result = st.session_state.get("query_expansion_result")
    if result:
        hpo_ids = result.get("hpo_ids", [])
        genes = result.get("genes", [])
        gene_count = result.get("gene_count", 0)
        
        st.success(f"Found {len(hpo_ids)} HPO terms and {gene_count} associated genes")
        
        # HPO Terms
        if hpo_ids:
            st.markdown("### Extracted HPO Terms")
            for hpo_id in hpo_ids:
                st.code(hpo_id, language=None)
        
        # Genes
        if genes:
            st.markdown(f"### Associated Genes ({len(genes)})")
            
            # Display in columns
            num_cols = 5
            cols = st.columns(num_cols)
            for i, gene in enumerate(genes[:50]):  # Limit to first 50 for display
                with cols[i % num_cols]:
                    st.markdown(f"- {gene}")
            
            if len(genes) > 50:
                st.caption(f"...and {len(genes) - 50} more genes")
            
            # Download
            df = pd.DataFrame({
                "Gene": genes,
                "HPO_Terms": [", ".join(hpo_ids)] * len(genes),
            })
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 Download Gene List (CSV)",
                data=csv,
                file_name="phenotype_genes.csv",
                mime="text/csv",
            )
        else:
            st.info("No genes found. Try adding HPO term IDs to your query (e.g., HP:0000001).")


if __name__ == "__main__":
    render_phenotypes_page()

