"""
Cross-omics dataset summary generation.

Generates multi-omics summaries for Datasets by aggregating evidence
from linked signatures, experiments, and programs.
"""

from __future__ import annotations

from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.context_extraction import (
    extract_aggregated_context,
    extract_disease_context,
    extract_matrix_context,
    extract_model_system_context,
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


def cross_omics_dataset_summary(
    dataset_page_id: str,
    top_k_chunks: int = 100,
) -> str:
    """
    Summarize cross-omics context for a single dataset.

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for dataset %s",
        dataset_page_id,
    )
    
    # Fetch dataset page
    dataset_page = fetch_notion_page(dataset_page_id)
    if not dataset_page:
        return f"Error: Could not fetch Dataset page {dataset_page_id}."
    
    dataset_name = extract_text_property(dataset_page, "Name") or "Unknown Dataset"
    omics_type = extract_select_values(dataset_page, "Omics Type")
    signature_ids = extract_relation_ids(dataset_page, "Related Signature(s)")
    experiment_ids = extract_relation_ids(dataset_page, "Related Experiments")
    program_ids = extract_relation_ids(dataset_page, "Related Programs")
    
    # Extract context from dataset page
    diseases = extract_disease_context(dataset_page)
    matrix = extract_matrix_context(dataset_page)
    model_systems = extract_model_system_context(dataset_page)
    
    dataset_context = {
        "diseases": diseases,
        "matrix": matrix,
        "model_systems": model_systems,
    }
    
    # Get dataset chunks
    dataset_chunks = retrieve_chunks_for_objects(
        [dataset_page_id],
        "dataset",
        top_k_per_object=top_k_chunks,
    )
    
    # Get signature chunks
    signature_chunks: List[Dict[str, Any]] = []
    for sig_id in signature_ids[:10]:  # Limit to 10 signatures
        sig_id_clean = sig_id.replace("-", "")
        chunks = query_pinecone(
            user_query="signature",
            top_k=5,
            meta_filter={"signature_page_id": sig_id_clean},
            source_types=["Signature"],
        )
        signature_chunks.extend(chunks)
    
    # Get experiment chunks
    experiment_chunks = retrieve_chunks_for_objects(
        experiment_ids,
        "experiment",
        top_k_per_object=10,
    )
    
    all_chunks = dataset_chunks + signature_chunks + experiment_chunks
    
    if not all_chunks:
        return (
            f"No sufficient multi-omics context found for dataset '{dataset_name}'. "
            "No chunks were found for this dataset or its related signatures/experiments."
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
    
    # Aggregate context from related pages
    all_related_ids = signature_ids + experiment_ids + program_ids
    aggregated_context = extract_aggregated_context(all_related_ids, page_type="dataset")
    
    # Merge dataset context with aggregated context
    if aggregated_context.get("diseases") or aggregated_context.get("matrix") or aggregated_context.get("model_systems"):
        dataset_context["diseases"] = list(set(dataset_context["diseases"] + aggregated_context.get("diseases", [])))
        dataset_context["matrix"] = list(set(dataset_context["matrix"] + aggregated_context.get("matrix", [])))
        dataset_context["model_systems"] = list(set(dataset_context["model_systems"] + aggregated_context.get("model_systems", [])))
    
    # Identify comparative context
    comparative_context = identify_comparative_context(dataset_context)
    
    # Build additional info
    additional_info = f"""Dataset details:
- Omics Type: {', '.join(omics_type) if omics_type else 'Not specified'}
- Linked signatures: {len(signature_ids)}
- Linked experiments: {len(experiment_ids)}
- Linked programs: {len(program_ids)}
- Signature chunks: {len(signature_chunks)}
- Experiment chunks: {len(experiment_chunks)}"""
    
    # Build enhanced prompt
    prompt = build_enhanced_prompt(
        entity_name=dataset_name,
        entity_type="dataset",
        context_info=dataset_context if (dataset_context.get("diseases") or dataset_context.get("matrix") or dataset_context.get("model_systems")) else None,
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

