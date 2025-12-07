"""
Postgres-based cross-omics feature summary generation.

Generates multi-omics summaries for Features (genes, proteins, metabolites, lipids)
by aggregating evidence from linked datasets stored in Postgres.

This is the Postgres-first version that:
- Queries Postgres directly (no Notion API calls required)
- Uses `notion_page_id` from Postgres models for chunk retrieval (temporary bridge
  until Pinecone is re-indexed with Postgres UUIDs)
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Feature as FeatureModel
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
from amprenta_rag.query.pinecone_query import query_pinecone

logger = get_logger(__name__)


def cross_omics_feature_summary_postgres(
    feature_id: UUID | None = None,
    feature_name: str | None = None,
    feature_type: str | None = None,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Generate a cross-omics summary for a Feature using Postgres data.

    Can find feature by UUID or by name+type. Queries Postgres for the feature
    and its linked datasets, then retrieves chunks from Pinecone.

    Args:
        feature_id: Postgres UUID of Feature (if provided, name/type ignored)
        feature_name: Feature name (e.g., "TP53", "Cer(d18:1/16:0)")
        feature_type: Feature type ("gene", "protein", "metabolite", "lipid")
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS][POSTGRES] Generating cross-omics summary for feature",
    )

    # Get database session
    db = next(get_db())
    try:
        # Find feature in Postgres
        feature: FeatureModel | None = None

        if feature_id:
            # Query by UUID
            feature = db.query(FeatureModel).filter(FeatureModel.id == feature_id).first()
            if not feature:
                return f"Error: Feature with ID {feature_id} not found in Postgres."
        elif feature_name and feature_type:
            # Query by name and type
            feature = (
                db.query(FeatureModel)
                .filter(
                    FeatureModel.name == feature_name,
                    FeatureModel.feature_type == feature_type.lower(),
                )
                .first()
            )
            if not feature:
                # Try normalized name
                feature = (
                    db.query(FeatureModel)
                    .filter(
                        FeatureModel.normalized_name == feature_name,
                        FeatureModel.feature_type == feature_type.lower(),
                    )
                    .first()
                )
            if not feature:
                return f"Error: Feature '{feature_name}' of type '{feature_type}' " "not found in Postgres."
        else:
            return "Error: Must provide either feature_id or (feature_name and feature_type)."

        feature_name_display = feature.name or "Unknown Feature"
        feature_type_display = feature.feature_type or "unknown"

        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found feature: %s (%s)",
            feature_name_display,
            feature_type_display,
        )

        # Get linked datasets via relationship
        datasets = feature.datasets

        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found %d linked datasets for feature %s",
            len(datasets),
            feature_name_display,
        )

        # Get Notion page IDs for chunk retrieval (backward compatibility)
        dataset_page_ids: List[str] = []
        for dataset in datasets[:top_k_datasets]:
            if dataset.notion_page_id:
                dataset_page_ids.append(dataset.notion_page_id)

        # Query Pinecone for chunks mentioning this feature
        query_text = f"{feature_type_display} {feature_name_display} multi-omics"
        feature_chunks = query_pinecone(
            user_query=query_text,
            top_k=top_k_chunks,
            meta_filter=None,
            source_types=None,
        )

        # Filter chunks that actually mention the feature
        filtered_chunks: List[Dict[str, Any]] = []
        feature_name_lower = feature_name_display.lower()
        for chunk in feature_chunks:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            snippet = meta.get("snippet", "").lower()
            title = meta.get("title", "").lower()
            text = f"{title} {snippet}"

            if feature_name_lower in text:
                filtered_chunks.append(chunk)

        # Also get chunks from linked datasets
        if dataset_page_ids:
            dataset_chunks = retrieve_chunks_for_objects(
                dataset_page_ids,
                "dataset",
                top_k_per_object=top_k_chunks // max(len(dataset_page_ids), 1),
            )
            filtered_chunks.extend(dataset_chunks)

        if not filtered_chunks:
            logger.warning(
                "[RAG][CROSS-OMICS][POSTGRES] No chunks found for feature %s. "
                "Linked datasets may not have been ingested into RAG.",
                feature_name_display,
            )

            # Fall back to metadata-only summary
            aggregated_context = aggregate_context_from_models(datasets, [])
            comparative_context = identify_comparative_context_postgres(aggregated_context)

            omics_types = {}
            for dataset in datasets:
                omics_type = dataset.omics_type or "Other"
                omics_types[omics_type] = omics_types.get(omics_type, 0) + 1

            additional_info = f"""This {feature_type_display} feature appears in:
- {len(datasets)} linked dataset(s)
- Omics types: {', '.join(omics_types.keys()) if omics_types else 'None'}

Summarize how this feature behaves across different omics modalities and contexts."""

            prompt = build_enhanced_prompt(
                entity_name=feature_name_display,
                entity_type=f"{feature_type_display} feature",
                context_info=(
                    aggregated_context
                    if aggregated_context.get("diseases")
                    or aggregated_context.get("matrix")
                    or aggregated_context.get("model_systems")
                    else None
                ),
                omics_counts=omics_types,
                additional_info=additional_info,
            )

            summary = synthesize_cross_omics_summary(
                prompt,
                [],
                include_comparative=comparative_context is not None,
            )

            return summary

        # Group by omics type
        chunks_by_omics = group_chunks_by_omics_type(filtered_chunks)

        # Get chunk texts
        context_chunks: List[str] = []
        for chunk in filtered_chunks[:top_k_chunks]:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            title = meta.get("title", "")
            source = meta.get("source_type", "")

            chunk_text = get_chunk_text(chunk)
            if chunk_text:
                context_chunks.append(f"[{source}] {title}\n{chunk_text}")

        omics_counts = {omics: len(chunks) for omics, chunks in chunks_by_omics.items() if chunks}

        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Retrieved %d chunks for %s feature '%s': %s",
            len(filtered_chunks),
            feature_type_display,
            feature_name_display,
            omics_counts,
        )

        # Aggregate context from linked datasets
        aggregated_context = aggregate_context_from_models(datasets, [])
        comparative_context = identify_comparative_context_postgres(aggregated_context)

        # Build additional info
        additional_info = f"""This {feature_type_display} feature appears in:
- {len(datasets)} linked dataset(s)

Summarize how this feature behaves across different omics modalities and contexts."""

        # Build enhanced prompt
        prompt = build_enhanced_prompt(
            entity_name=feature_name_display,
            entity_type=f"{feature_type_display} feature",
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

    finally:
        db.close()
