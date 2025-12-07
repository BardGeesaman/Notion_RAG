"""
Postgres-based cross-omics program summary generation.

Generates multi-omics summaries for Programs by aggregating evidence
from linked experiments and datasets stored in Postgres.

This is the Postgres-first version that:
- Queries Postgres directly (no Notion API calls required)
- Uses `notion_page_id` from Postgres models for chunk retrieval (temporary bridge
  until Pinecone is re-indexed with Postgres UUIDs)
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.database.models import Experiment as ExperimentModel
from amprenta_rag.database.models import Program as ProgramModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.helpers import (
    get_chunk_text,
    group_chunks_by_omics_type,
    retrieve_chunks_for_objects,
)
from amprenta_rag.query.cross_omics.prompt_templates import build_enhanced_prompt
from amprenta_rag.query.cross_omics.synthesis import synthesize_cross_omics_summary

logger = get_logger(__name__)


def extract_context_from_postgres_model(
    model: DatasetModel | ExperimentModel,
) -> Dict[str, Any]:
    """
    Extract context information from a Postgres model.

    Args:
        model: Dataset or Experiment model instance

    Returns:
        Dictionary with disease, matrix, model_systems
    """
    context: Dict[str, Any] = {
        "diseases": model.disease or [],
        "matrix": [],
        "model_systems": [],
    }

    # Extract matrix from Dataset or Experiment
    if hasattr(model, "matrix") and model.matrix:
        context["matrix"] = model.matrix

    # Extract model systems from Experiment
    if hasattr(model, "model_systems") and model.model_systems:
        context["model_systems"] = model.model_systems

    return context


def aggregate_context_from_models(
    datasets: List[DatasetModel],
    experiments: List[ExperimentModel],
) -> Dict[str, Any]:
    """
    Aggregate context from multiple Postgres models.

    Args:
        datasets: List of Dataset models
        experiments: List of Experiment models

    Returns:
        Aggregated context dictionary
    """
    all_diseases: List[str] = []
    all_matrix: List[str] = []
    all_model_systems: List[str] = []

    for dataset in datasets:
        if dataset.disease:
            all_diseases.extend(dataset.disease)
        if hasattr(dataset, "matrix") and dataset.matrix:
            all_matrix.extend(dataset.matrix)

    for experiment in experiments:
        if experiment.disease:
            all_diseases.extend(experiment.disease)
        if experiment.matrix:
            all_matrix.extend(experiment.matrix)
        if experiment.model_systems:
            all_model_systems.extend(experiment.model_systems)

    return {
        "diseases": list(set(all_diseases)),
        "matrix": list(set(all_matrix)),
        "model_systems": list(set(all_model_systems)),
    }


def identify_comparative_context_postgres(
    aggregated_context: Dict[str, Any],
) -> Dict[str, Any] | None:
    """
    Identify if there's comparative context (multiple diseases, matrices, etc.).

    Args:
        aggregated_context: Aggregated context dictionary

    Returns:
        Comparative context dict if applicable, None otherwise
    """
    diseases = aggregated_context.get("diseases", [])
    matrices = aggregated_context.get("matrix", [])
    model_systems = aggregated_context.get("model_systems", [])

    if len(diseases) > 1 or len(matrices) > 1 or len(model_systems) > 1:
        return {
            "multiple_diseases": len(diseases) > 1,
            "multiple_matrices": len(matrices) > 1,
            "multiple_model_systems": len(model_systems) > 1,
            "diseases": diseases,
            "matrices": matrices,
            "model_systems": model_systems,
        }

    return None


def cross_omics_program_summary_postgres(
    program_id: UUID,
    top_k_per_omics: int = 20,
) -> str:
    """
    Generate a cross-omics summary for a Program using Postgres data.

    This function queries Postgres for the program and its linked
    experiments/datasets, then retrieves chunks from Pinecone.

    Args:
        program_id: Postgres UUID of Program
        top_k_per_omics: Maximum chunks to retrieve per omics type

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS][POSTGRES] Generating cross-omics summary for program %s",
        program_id,
    )

    # Get database session
    db = next(get_db())
    try:
        # Query program from Postgres
        program = db.query(ProgramModel).filter(ProgramModel.id == program_id).first()
        if not program:
            return f"Error: Program with ID {program_id} not found in Postgres."

        program_name = program.name or "Unknown Program"

        # Get linked experiments and datasets via relationships
        experiments = program.experiments
        datasets = program.datasets

        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Found %d experiments, %d datasets for program %s",
            len(experiments),
            len(datasets),
            program_name,
        )

        if not experiments and not datasets:
            return (
                f"No sufficient multi-omics context found for program '{program_name}'. "
                "No experiments or datasets are linked to this program."
            )

        # Get Notion page IDs for chunk retrieval (backward compatibility)
        # Chunks in Pinecone are currently indexed with Notion page IDs
        dataset_page_ids: List[str] = []
        experiment_page_ids: List[str] = []

        for dataset in datasets:
            if dataset.notion_page_id:
                dataset_page_ids.append(dataset.notion_page_id)

        for experiment in experiments:
            if experiment.notion_page_id:
                experiment_page_ids.append(experiment.notion_page_id)

        # If no Notion page IDs, we can't retrieve chunks yet
        # Note: Chunks in Pinecone are currently indexed with Notion page IDs.
        # Once Pinecone is re-indexed with Postgres UUIDs, we can retrieve chunks
        # directly using Postgres UUIDs instead of notion_page_id.
        if not dataset_page_ids and not experiment_page_ids:
            logger.warning(
                "[RAG][CROSS-OMICS][POSTGRES] No Notion page IDs found for datasets/experiments. "
                "Chunks may not be retrievable. Consider re-indexing with Postgres UUIDs."
            )
            # Still generate summary from metadata
            aggregated_context = aggregate_context_from_models(datasets, experiments)
            comparative_context = identify_comparative_context_postgres(aggregated_context)

            # Build summary from metadata only
            omics_types = {}
            for dataset in datasets:
                omics_type = dataset.omics_type or "Other"
                omics_types[omics_type] = omics_types.get(omics_type, 0) + 1

            additional_info = f"""The program is linked to:
- {len(experiments)} experiment(s)
- {len(datasets)} dataset(s)
- Omics types: {', '.join(omics_types.keys()) if omics_types else 'None'}"""

            prompt = build_enhanced_prompt(
                entity_name=program_name,
                entity_type="program",
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

            # Generate summary from metadata only (no chunks)
            summary = synthesize_cross_omics_summary(
                prompt,
                [],  # No chunks available
                include_comparative=comparative_context is not None,
            )

            return summary

        # Retrieve chunks from datasets
        dataset_chunks = retrieve_chunks_for_objects(
            dataset_page_ids,
            "dataset",
            top_k_per_object=top_k_per_omics,
        )

        # Retrieve chunks from experiments
        experiment_chunks = retrieve_chunks_for_objects(
            experiment_page_ids,
            "experiment",
            top_k_per_object=top_k_per_omics,
        )

        all_chunks = dataset_chunks + experiment_chunks

        if not all_chunks:
            logger.warning(
                "[RAG][CROSS-OMICS][POSTGRES] No chunks found for program %s. "
                "Linked datasets and experiments may not have been ingested into RAG.",
                program_name,
            )
            # Fall back to metadata-only summary
            aggregated_context = aggregate_context_from_models(datasets, experiments)
            comparative_context = identify_comparative_context_postgres(aggregated_context)

            omics_types = {}
            for dataset in datasets:
                omics_type = dataset.omics_type or "Other"
                omics_types[omics_type] = omics_types.get(omics_type, 0) + 1

            additional_info = f"""The program is linked to:
- {len(experiments)} experiment(s)
- {len(datasets)} dataset(s)
- Omics types: {', '.join(omics_types.keys()) if omics_types else 'None'}"""

            prompt = build_enhanced_prompt(
                entity_name=program_name,
                entity_type="program",
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

        # Extract aggregated context from Postgres models
        aggregated_context = aggregate_context_from_models(datasets, experiments)
        comparative_context = identify_comparative_context_postgres(aggregated_context)

        # Group chunks by omics type
        chunks_by_omics = group_chunks_by_omics_type(all_chunks)

        # Get chunk texts
        context_chunks: List[str] = []
        for chunk in all_chunks:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            title = meta.get("title", "")
            source = meta.get("source_type", "")

            chunk_text = get_chunk_text(chunk)
            if chunk_text:
                context_chunks.append(f"[{source}] {title}\n{chunk_text}")

        omics_counts = {omics: len(chunks) for omics, chunks in chunks_by_omics.items() if chunks}

        logger.info(
            "[RAG][CROSS-OMICS][POSTGRES] Retrieved %d total chunks for program %s: %s",
            len(all_chunks),
            program_name,
            omics_counts,
        )

        # Build additional info
        additional_info = f"""The program is linked to:
- {len(experiments)} experiment(s)
- {len(datasets)} dataset(s)"""

        # Build enhanced prompt with context
        prompt = build_enhanced_prompt(
            entity_name=program_name,
            entity_type="program",
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
