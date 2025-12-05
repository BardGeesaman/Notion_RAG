"""
Cross-omics signature summary generation.

Generates multi-omics summaries for Signatures by aggregating evidence
from matched datasets and signature components.
"""

from __future__ import annotations

from typing import Any, Dict, List

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.context_extraction import (
    extract_aggregated_context,
    extract_disease_context,
    extract_matrix_context,
    identify_comparative_context,
)
from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_select_values,
    extract_text_property,
    fetch_notion_page,
    get_chunk_text,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt
from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary
from amprenta_rag.query.pinecone_query import query_pinecone

logger = get_logger(__name__)


def cross_omics_signature_summary(
    signature_page_id: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Generate a cross-omics summary for a Signature.

    Args:
        signature_page_id: Notion page ID of Signature (with dashes)
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for signature %s",
        signature_page_id,
    )
    
    # Fetch signature page
    signature_page = fetch_notion_page(signature_page_id)
    if not signature_page:
        return f"Error: Could not fetch Signature page {signature_page_id}."
    
    signature_name = extract_text_property(signature_page, "Name") or "Unknown Signature"
    modalities = extract_select_values(signature_page, "Modalities")
    
    # Extract context from signature page
    diseases = extract_disease_context(signature_page)
    matrix = extract_matrix_context(signature_page)
    
    signature_context = {
        "diseases": diseases,
        "matrix": matrix,
        "model_systems": [],  # Signatures typically don't have model systems
    }
    
    # Find datasets with this signature
    # Query Experimental Data Assets database for datasets with this signature in Related Signature(s)
    cfg = get_config()
    dataset_ids: List[str] = []
    
    try:
        # Query Pinecone for chunks with this signature
        signature_id_clean = signature_page_id.replace("-", "")
        meta_filter = {"signature_page_id": signature_id_clean}
        
        signature_chunks = query_pinecone(
            user_query=f"signature {signature_name} multi-omics analysis",
            top_k=top_k_chunks,
            meta_filter=meta_filter,
            source_types=["Signature", "Dataset"],
        )
        
        # Extract unique dataset IDs from chunks
        for chunk in signature_chunks:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            dataset_id = meta.get("dataset_page_id")
            if dataset_id:
                # Add dashes back
                dataset_id_with_dashes = f"{dataset_id[:8]}-{dataset_id[8:12]}-{dataset_id[12:16]}-{dataset_id[16:20]}-{dataset_id[20:]}"
                if dataset_id_with_dashes not in dataset_ids:
                    dataset_ids.append(dataset_id_with_dashes)
        
        # Also query Experimental Data Assets database directly
        try:
            exp_data_db_id = cfg.notion.experimental_data_assets_db_id
            if exp_data_db_id:
                url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
                payload = {
                    "filter": {
                        "property": "Related Signature(s)",
                        "relation": {"contains": signature_page_id},
                    },
                    "page_size": 100,
                }
                
                resp = requests.post(
                    url,
                    headers=notion_headers(),
                    json=payload,
                    timeout=30,
                )
                resp.raise_for_status()
                
                data = resp.json()
                db_dataset_ids = [r.get("id", "") for r in data.get("results", []) if r.get("id")]
                dataset_ids.extend(db_dataset_ids)
                dataset_ids = list(set(dataset_ids))  # Deduplicate
                
        except Exception as e:
            logger.debug(
                "[RAG][CROSS-OMICS] Could not query Experimental Data Assets DB: %r",
                e,
            )
    
    except Exception as e:
        logger.warning(
            "[RAG][CROSS-OMICS] Error finding datasets for signature: %r",
            e,
        )
    
    logger.info(
        "[RAG][CROSS-OMICS] Found %d datasets for signature %s",
        len(dataset_ids),
        signature_name,
    )
    
    # Limit to top_k_datasets
    dataset_ids = dataset_ids[:top_k_datasets]
    
    # Retrieve chunks from datasets and signature
    all_chunks: List[Dict[str, Any]] = []
    
    # Get signature chunks
    signature_id_clean = signature_page_id.replace("-", "")
    signature_chunks = query_pinecone(
        user_query=f"signature {signature_name}",
        top_k=20,
        meta_filter={"signature_page_id": signature_id_clean},
        source_types=["Signature"],
    )
    all_chunks.extend(signature_chunks)
    
    # Get dataset chunks
    if dataset_ids:
        dataset_chunks = retrieve_chunks_for_objects(
            dataset_ids,
            "dataset",
            top_k_per_object=top_k_chunks // max(len(dataset_ids), 1),
        )
        all_chunks.extend(dataset_chunks)
    
    if not all_chunks:
        return (
            f"No sufficient multi-omics context found for signature '{signature_name}'. "
            "No chunks were found for this signature or its matched datasets."
        )
    
    # Group by omics type
    chunks_by_omics = group_chunks_by_omics_type(all_chunks)
    
    # Get chunk texts
    context_chunks: List[str] = []
    for chunk in all_chunks[:top_k_chunks]:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        chunk_text = get_chunk_text(chunk)
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    # Aggregate context from matched datasets
    if dataset_ids:
        aggregated_context = extract_aggregated_context(dataset_ids, page_type="dataset")
        # Merge signature context with aggregated context
        if aggregated_context.get("diseases") or aggregated_context.get("matrix") or aggregated_context.get("model_systems"):
            signature_context["diseases"] = list(set(signature_context["diseases"] + aggregated_context.get("diseases", [])))
            signature_context["matrix"] = list(set(signature_context["matrix"] + aggregated_context.get("matrix", [])))
            signature_context["model_systems"] = list(set(signature_context["model_systems"] + aggregated_context.get("model_systems", [])))
    
    # Identify comparative context
    comparative_context = identify_comparative_context(signature_context)
    
    # Build additional info
    additional_info = f"""Signature details:
- Modalities: {', '.join(modalities) if modalities else 'Not specified'}
- Matched datasets: {len(dataset_ids)}
- Signature chunks: {len(signature_chunks)}"""
    
    # Build enhanced prompt
    prompt = build_enhanced_prompt(
        entity_name=signature_name,
        entity_type="signature",
        context_info=signature_context if (signature_context.get("diseases") or signature_context.get("matrix")) else None,
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

