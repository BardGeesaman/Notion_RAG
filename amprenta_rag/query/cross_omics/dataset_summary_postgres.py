"""
Postgres-based cross-omics dataset summary generation.

Generates multi-omics summaries for Datasets by aggregating evidence
from the dataset itself and related programs/experiments stored in Postgres.

This is the Postgres-first version that:
- Queries Postgres directly (no Notion API calls required)
- Uses `notion_page_id` from Postgres models for chunk retrieval (temporary bridge
  until Pinecone is re-indexed with Postgres UUIDs)
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.helpers import (
    get_chunk_text,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
from amprenta_rag.query.cross_omics.program_summary_postgres import (
    aggregate_context_from_models,
    identify_comparative_context_postgres,
)
from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt
from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary

logger = get_logger(__name__)


def cross_omics_dataset_summary_postgres(
    dataset_id: UUID,
    top_k_chunks: int = 20,
) -> str:
    """
    Generate a cross-omics summary for a Dataset using Postgres data.

    Queries Postgres for the dataset and its metadata, then retrieves
    chunks from Pinecone.

    Args:
        dataset_id: Postgres UUID of Dataset
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS][POSTGRES] Generating cross-omics summary for dataset %s",
        dataset_id,
    )

    with db_session() as db:
        # Query dataset from Postgres
        dataset = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
        if not dataset:
            return f"Error: Dataset with ID {dataset_id} not found in Postgres."
        
        dataset_name = dataset.name or "Unknown Dataset"
        omics_type = dataset.omics_type or "Unknown"
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found dataset: %s (%s)",
            dataset_name,
            omics_type,
        )
        
        programs = dataset.programs
        experiments = dataset.experiments
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Dataset linked to %d program(s), %d experiment(s)",
            len(programs),
            len(experiments),
        )
        
        dataset_page_ids: List[str] = []
        if dataset.notion_page_id:
            dataset_page_ids.append(dataset.notion_page_id)
        
        all_chunks: List[Dict[str, Any]] = []
        
        if dataset_page_ids:
            dataset_chunks = retrieve_chunks_for_objects(
                dataset_page_ids,
                "dataset",
                top_k_per_object=top_k_chunks,
            )
            all_chunks.extend(dataset_chunks)
        
        if not all_chunks:
            logger.warning(
                "[RAG][CROSS-OMICS][POSTGRES] No chunks found for dataset %s. "
                "Dataset may not have been ingested into RAG.",
                dataset_name,
            )
            
            aggregated_context = aggregate_context_from_models([dataset], experiments)
            comparative_context = identify_comparative_context_postgres(aggregated_context)
            
            omics_counts = {omics_type: 1} if omics_type != "Unknown" else {}
            
            additional_info = f"""This dataset:
- Omics type: {omics_type}
- Linked to {len(programs)} program(s)
- Linked to {len(experiments)} experiment(s)
- Description: {dataset.description or 'No description available'}"""

            prompt = build_enhanced_prompt(
                entity_name=dataset_name,
                entity_type="dataset",
                context_info=(
                    aggregated_context
                    if aggregated_context.get("diseases")
                    or aggregated_context.get("matrix")
                    or aggregated_context.get("model_systems")
                    else None
                ),
                omics_counts=omics_counts,
                additional_info=additional_info,
            )

            summary = synthesize_cross_omics_summary(
                prompt,
                [],
                include_comparative=comparative_context is not None,
            )

            return summary

        chunks_by_omics = group_chunks_by_omics_type(all_chunks)
        
        context_chunks: List[str] = []
        for chunk in all_chunks[:top_k_chunks]:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            title = meta.get("title", "")
            source = meta.get("source_type", "")
            
            chunk_text = get_chunk_text(chunk)
            if chunk_text:
                context_chunks.append(f"[{source}] {title}\n{chunk_text}")
        
        omics_counts = {omics: len(chunks) for omics, chunks in chunks_by_omics.items() if chunks}
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Retrieved %d chunks for dataset %s: %s",
            len(all_chunks),
            dataset_name,
            omics_counts,
        )
        
        aggregated_context = aggregate_context_from_models([dataset], experiments)
        comparative_context = identify_comparative_context_postgres(aggregated_context)
        
        additional_info = f"""This dataset:
- Omics type: {omics_type}
- Linked to {len(programs)} program(s)
- Linked to {len(experiments)} experiment(s)
- Description: {dataset.description or 'No description available'}"""

        # Build enhanced prompt
        prompt = build_enhanced_prompt(
            entity_name=dataset_name,
            entity_type="dataset",
            context_info=(
                aggregated_context
                if aggregated_context
                and (
                    aggregated_context.get("diseases")
                    or aggregated_context.get("matrix")
                    or aggregated_context.get("model_systems")
                )
                else None
            ),
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
