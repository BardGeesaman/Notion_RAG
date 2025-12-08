"""
Cross-omics feature summary generation (Postgres-backed wrapper).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import get_page_text, notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[CROSS-OMICS][FEATURE-SUMMARY] notion_headers() deprecated - Notion support removed")
    return {}

def get_page_text(page_id: str) -> str:
    """DEPRECATED: Notion support removed. Returns empty string."""
    logger.debug("[CROSS-OMICS][FEATURE-SUMMARY] get_page_text() deprecated - Notion support removed")
    return ""
from amprenta_rag.query.cross_omics.context_extraction import (
    extract_aggregated_context,
    identify_comparative_context,
)
from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    fetch_notion_page,
    get_chunk_text,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt
from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary
from amprenta_rag.query.pinecone_query import query_pinecone

logger = get_logger(__name__)


def cross_omics_feature_summary(
    feature_name: str,
    feature_type: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Summarize all multi-omics evidence for a single feature.

    Args:
        feature_name: Feature name (e.g., "TP53", "Cer(d18:1/16:0)")
        feature_type: "gene", "protein", "metabolite", or "lipid"
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for %s feature: %s",
        feature_type,
        feature_name,
    )
    
    # Find feature page in appropriate database
    cfg = get_config()
    
    # Map feature type to database ID
    db_map = {
        "gene": cfg.notion.gene_features_db_id if hasattr(cfg.notion, "gene_features_db_id") else None,
        "protein": cfg.notion.protein_features_db_id if hasattr(cfg.notion, "protein_features_db_id") else None,
        "metabolite": cfg.notion.metabolite_features_db_id if hasattr(cfg.notion, "metabolite_features_db_id") else None,
        "lipid": cfg.notion.lipid_species_db_id if hasattr(cfg.notion, "lipid_species_db_id") else None,
    }
    
    db_id = db_map.get(feature_type)
    if not db_id:
        return f"Error: {feature_type} Features database ID not configured."
    
    # Query for feature page
    feature_page_id: Optional[str] = None
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Name",
                "title": {"equals": feature_name},
            },
            "page_size": 1,
        }
        
        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        
        results = resp.json().get("results", [])
        if results:
            feature_page_id = results[0].get("id", "")
    except Exception as e:
        logger.warning(
            "[RAG][CROSS-OMICS] Error finding %s feature page for '%s': %r",
            feature_type,
            feature_name,
            e,
        )
    
    # Find linked datasets, experiments, programs
    dataset_ids: List[str] = []
    
    if feature_page_id:
        feature_page = fetch_notion_page(feature_page_id)
        if feature_page:
            # Extract relations based on feature type
            relation_property_candidates = {
                "gene": ["Transcriptomics Datasets", "Datasets", "Related Datasets"],
                "protein": ["Proteomics Datasets", "Datasets", "Related Datasets"],
                "metabolite": ["Metabolomics Datasets", "Datasets", "Related Datasets"],
                "lipid": ["Experimental Data Assets", "Datasets", "Related Datasets"],
            }
            
            candidates = relation_property_candidates.get(feature_type, ["Datasets", "Related Datasets"])
            for prop_name in candidates:
                dataset_ids = extract_relation_ids(feature_page, prop_name)
                if dataset_ids:
                    break
    
    # Query Pinecone for chunks mentioning this feature
    query_text = f"{feature_type} {feature_name} multi-omics"
    feature_chunks = query_pinecone(
        user_query=query_text,
        top_k=top_k_chunks,
        meta_filter=None,
        source_types=None,
    )
    
    # Filter chunks that actually mention the feature
    filtered_chunks = []
    feature_name_lower = feature_name.lower()
    for chunk in feature_chunks:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        snippet = meta.get("snippet", "").lower()
        title = meta.get("title", "").lower()
        text = f"{title} {snippet}"
        
        if feature_name_lower in text:
            filtered_chunks.append(chunk)
    
    # Also get chunks from linked datasets
    if dataset_ids:
        dataset_chunks = retrieve_chunks_for_objects(
            dataset_ids[:top_k_datasets],
            "dataset",
            top_k_per_object=top_k_chunks // max(len(dataset_ids[:top_k_datasets]), 1),
        )
        filtered_chunks.extend(dataset_chunks)
    
    if not filtered_chunks:
        return (
            f"No sufficient multi-omics context found for {feature_type} feature '{feature_name}'. "
            "No relevant chunks or linked datasets were found."
        )
    
    # Group by omics type
    chunks_by_omics = group_chunks_by_omics_type(filtered_chunks)
    
    # Get chunk texts - try to get full text from Notion, fallback to snippet
    context_chunks: List[str] = []
    for chunk in filtered_chunks[:top_k_chunks]:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        # Try to get full chunk text from Notion
        chunk_text = ""
        chunk_page_id = meta.get("notion_chunk_page_id")
        if chunk_page_id:
            try:
                full_text = get_page_text(chunk_page_id)
                if full_text:
                    chunk_text = full_text
            except Exception:
                pass
        
        # Fallback to snippet if no full text
        if not chunk_text:
            snippet = meta.get("snippet", "")
            chunk_text = snippet
        
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    logger.info(
        "[RAG][CROSS-OMICS] Retrieved %d chunks for %s feature '%s': %s",
        len(filtered_chunks),
        feature_type,
        feature_name,
        omics_counts,
    )
    
    # Aggregate context from linked datasets
    aggregated_context = None
    if dataset_ids:
        aggregated_context = extract_aggregated_context(dataset_ids, page_type="dataset")
    
    # Identify comparative context
    comparative_context = None
    if aggregated_context:
        comparative_context = identify_comparative_context(aggregated_context)
    
    # Build additional info
    additional_info = f"""This {feature_type} feature appears in:
- {len(dataset_ids)} linked dataset(s)

Summarize how this feature behaves across different omics modalities and contexts."""
    
    # Build enhanced prompt
    prompt = build_enhanced_prompt(
        entity_name=feature_name,
        entity_type=f"{feature_type} feature",
        context_info=aggregated_context if aggregated_context and (aggregated_context.get("diseases") or aggregated_context.get("matrix") or aggregated_context.get("model_systems")) else None,
        omics_counts=omics_counts,
        additional_info=additional_info,
    )
    
    # Synthesize summary with comparative analysis if applicable
    include_comparative = comparative_context is not None
    summary = synthesize_cross_omics_summary(
        prompt,
        context_chunks,
        include_comparative=include_comparative,
    )
    
    return summary

