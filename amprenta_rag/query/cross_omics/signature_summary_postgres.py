"""
Postgres-based cross-omics signature summary generation.

Generates multi-omics summaries for Signatures by aggregating evidence
from linked datasets and signature components stored in Postgres.

This is the Postgres-first version that:
- Queries Postgres directly (no Notion API calls required)
- Uses `notion_page_id` from Postgres models for chunk retrieval (temporary bridge
  until Pinecone is re-indexed with Postgres UUIDs)
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Signature as SignatureModel
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


def cross_omics_signature_summary_postgres(
    signature_id: UUID | None = None,
    signature_name: str | None = None,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Generate a cross-omics summary for a Signature using Postgres data.

    Can find signature by UUID or by name. Queries Postgres for the signature
    and its linked datasets/components, then retrieves chunks from Pinecone.

    Args:
        signature_id: Postgres UUID of Signature (if provided, name ignored)
        signature_name: Signature name (e.g., "ALS-CSF-Core-6Ceramides")
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS][POSTGRES] Generating cross-omics summary for signature",
    )

    with db_session() as db:
        signature: SignatureModel | None = None
        
        if signature_id:
            signature = db.query(SignatureModel).filter(SignatureModel.id == signature_id).first()
            if not signature:
                return f"Error: Signature with ID {signature_id} not found in Postgres."
        elif signature_name:
            signature = db.query(SignatureModel).filter(SignatureModel.name == signature_name).first()
            if not signature:
                signature = db.query(SignatureModel).filter(SignatureModel.short_id == signature_name).first()
            if not signature:
                return f"Error: Signature '{signature_name}' not found in Postgres."
        else:
            return "Error: Must provide either signature_id or signature_name."
        
        signature_name_display = signature.name or "Unknown Signature"
        modalities = signature.modalities or []
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found signature: %s (modalities: %s)",
            signature_name_display,
            modalities,
        )
        
        datasets = signature.datasets
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found %d linked datasets for signature %s",
            len(datasets),
            signature_name_display,
        )
        
        components = signature.components
        feature_types_in_signature = list(set([c.feature_type for c in components if c.feature_type]))
        
        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Signature has %d components across %s",
            len(components),
            feature_types_in_signature,
        )
        
        dataset_page_ids: List[str] = []
        for dataset in datasets[:top_k_datasets]:
            if dataset.notion_page_id:
                dataset_page_ids.append(dataset.notion_page_id)
        
        all_chunks: List[Dict[str, Any]] = []
        
        if signature.notion_page_id:
            signature_id_clean = signature.notion_page_id.replace("-", "")
            meta_filter = {"signature_page_id": signature_id_clean}
            
            signature_chunks = query_pinecone(
                user_query=f"signature {signature_name_display} multi-omics analysis",
                top_k=top_k_chunks,
                meta_filter=meta_filter,
                source_types=["Signature", "Dataset"],
            )
            all_chunks.extend(signature_chunks)
        
        if dataset_page_ids:
            dataset_chunks = retrieve_chunks_for_objects(
                dataset_page_ids,
                "dataset",
                top_k_per_object=top_k_chunks // max(len(dataset_page_ids), 1),
            )
            all_chunks.extend(dataset_chunks)
        
        if not all_chunks:
            logger.warning(
                "[RAG][CROSS-OMICS][POSTGRES] No chunks found for signature %s. "
                "Linked datasets may not have been ingested into RAG.",
                signature_name_display,
            )
            
            aggregated_context = aggregate_context_from_models(datasets, [])
            comparative_context = identify_comparative_context_postgres(aggregated_context)
            
            omics_types = {}
            for dataset in datasets:
                omics_type = dataset.omics_type or "Other"
                omics_types[omics_type] = omics_types.get(omics_type, 0) + 1
            
            component_summary = []
            for comp in components[:10]:
                comp_str = comp.feature_name or f"{comp.feature_type} feature"
                if comp.direction:
                    comp_str += f" ({comp.direction})"
                if comp.weight:
                    comp_str += f" [weight: {comp.weight}]"
                component_summary.append(comp_str)
            
            additional_info = f"""This signature contains:
- {len(components)} component(s) across {len(feature_types_in_signature)} feature type(s): {', '.join(feature_types_in_signature)}
- {len(datasets)} linked dataset(s)
- Omics types: {', '.join(omics_types.keys()) if omics_types else 'None'}

Key components: {', '.join(component_summary[:5])}{'...' if len(component_summary) > 5 else ''}"""

            prompt = build_enhanced_prompt(
                entity_name=signature_name_display,
                entity_type="signature",
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
            "[RAG][CROSS-OMICS][POSTGRES] Retrieved %d chunks for signature %s: %s",
            len(all_chunks),
            signature_name_display,
            omics_counts,
        )
        
        aggregated_context = aggregate_context_from_models(datasets, [])
        comparative_context = identify_comparative_context_postgres(aggregated_context)
        
        component_summary = []
        for comp in components[:10]:
            comp_str = comp.feature_name or f"{comp.feature_type} feature"
            if comp.direction:
                comp_str += f" ({comp.direction})"
            if comp.weight:
                comp_str += f" [weight: {comp.weight}]"
            component_summary.append(comp_str)
        
        additional_info = f"""This signature contains:
- {len(components)} component(s) across {len(feature_types_in_signature)} feature type(s): {', '.join(feature_types_in_signature)}
- {len(datasets)} linked dataset(s)

Key components: {', '.join(component_summary[:5])}{'...' if len(component_summary) > 5 else ''}"""

        # Build enhanced prompt
        prompt = build_enhanced_prompt(
            entity_name=signature_name_display,
            entity_type="signature",
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
