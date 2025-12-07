"""
Cross-omics program summary generation.

Generates multi-omics summaries for Programs by aggregating evidence
from linked experiments and datasets.
"""

from __future__ import annotations

from typing import Any, Dict, List

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import get_page_text

def get_page_text(page_id: str) -> str:
    """DEPRECATED: Notion support removed. Returns empty string."""
    logger.debug("[CROSS-OMICS][PROGRAM-SUMMARY] get_page_text() deprecated - Notion support removed")
    return ""
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.context_extraction import (
    extract_aggregated_context,
    identify_comparative_context,
)
from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_text_property,
    fetch_notion_page,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt
from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary

logger = get_logger(__name__)


def cross_omics_program_summary(
    program_page_id: str,
    top_k_per_omics: int = 20,
) -> str:
    """
    Generate a cross-omics summary for a Program.

    Args:
        program_page_id: Notion page ID of Program (with dashes)
        top_k_per_omics: Maximum chunks to retrieve per omics type

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for program %s",
        program_page_id,
    )
    
    # Fetch program page
    program_page = fetch_notion_page(program_page_id)
    if not program_page:
        return f"Error: Could not fetch Program page {program_page_id}."
    
    # Try "Program" property first (actual schema), fallback to "Name"
    program_name = extract_text_property(program_page, "Program") or extract_text_property(program_page, "Name") or "Unknown Program"
    
    # Find related Experiments (schema shows "Experiments" relation)
    # Try "Experiments" first, then fallback to "Related Experiments"
    experiment_ids = extract_relation_ids(program_page, "Experiments")
    if not experiment_ids:
        experiment_ids = extract_relation_ids(program_page, "Related Experiments")
    
    # Datasets might be linked via Experiments, or check for "Program Datasets"
    dataset_ids = extract_relation_ids(program_page, "Program Datasets")
    
    logger.info(
        "[RAG][CROSS-OMICS] Found %d experiments, %d datasets for program %s",
        len(experiment_ids),
        len(dataset_ids),
        program_name,
    )
    
    if not experiment_ids and not dataset_ids:
        return (
            f"No sufficient multi-omics context found for program '{program_name}'. "
            "No experiments or datasets are linked to this program."
        )
    
    # Retrieve chunks from datasets
    dataset_chunks = retrieve_chunks_for_objects(
        dataset_ids,
        "dataset",
        top_k_per_object=top_k_per_omics,
    )
    
    # Retrieve chunks from experiments
    experiment_chunks = retrieve_chunks_for_objects(
        experiment_ids,
        "experiment",
        top_k_per_object=top_k_per_omics,
    )
    
    all_chunks = dataset_chunks + experiment_chunks
    
    if not all_chunks:
        return (
            f"No chunks found for program '{program_name}'. "
            "Linked datasets and experiments may not have been ingested into RAG."
        )
    
    # Extract aggregated context from all linked pages
    all_page_ids = dataset_ids + experiment_ids
    aggregated_context = extract_aggregated_context(all_page_ids, page_type="dataset")
    comparative_context = identify_comparative_context(aggregated_context)
    
    # Group chunks by omics type
    chunks_by_omics = group_chunks_by_omics_type(all_chunks)
    
    # Get chunk texts - try to get full text from Notion, fallback to snippet
    context_chunks: List[str] = []
    for chunk in all_chunks:
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
            except Exception as e:
                logger.debug(
                    "[RAG][CROSS-OMICS] Could not fetch full text for chunk %s: %r",
                    chunk_page_id,
                    e,
                )
        
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
        "[RAG][CROSS-OMICS] Retrieved %d total chunks for program %s: %s",
        len(all_chunks),
        program_name,
        omics_counts,
    )
    
    # Build additional info
    additional_info = f"""The program is linked to:
- {len(experiment_ids)} experiment(s)
- {len(dataset_ids)} dataset(s)"""
    
    # Build enhanced prompt with context
    prompt = build_enhanced_prompt(
        entity_name=program_name,
        entity_type="program",
        context_info=aggregated_context if aggregated_context.get("diseases") or aggregated_context.get("matrix") or aggregated_context.get("model_systems") else None,
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

