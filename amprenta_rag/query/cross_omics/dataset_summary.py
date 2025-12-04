"""
Cross-omics dataset summary generation.

Generates multi-omics summaries for Datasets by aggregating evidence
from linked signatures, experiments, and programs.
"""

from __future__ import annotations

from typing import Any, Dict, List

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_select_values,
    extract_text_property,
    fetch_notion_page,
    get_chunk_text,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
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
    
    # Build prompt
    prompt = f"""Generate a cross-omics summary for the dataset: {dataset_name}

Dataset details:
- Omics Type: {', '.join(omics_type) if omics_type else 'Not specified'}
- Linked signatures: {len(signature_ids)}
- Linked experiments: {len(experiment_ids)}
- Linked programs: {len(program_ids)}

Context chunks retrieved:
- {omics_counts.get('Lipidomics', 0)} lipidomics chunks
- {omics_counts.get('Metabolomics', 0)} metabolomics chunks
- {omics_counts.get('Proteomics', 0)} proteomics chunks
- {omics_counts.get('Transcriptomics', 0)} transcriptomics chunks
- {len(signature_chunks)} signature chunks
- {len(experiment_chunks)} experiment chunks
"""
    
    # Synthesize summary
    summary = synthesize_cross_omics_summary(prompt, context_chunks)
    
    return summary

