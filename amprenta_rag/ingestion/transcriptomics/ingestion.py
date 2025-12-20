"""
Main transcriptomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
transcriptomics dataset ingestion pipeline: file parsing, Postgres storage,
feature linking, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, cast
from uuid import UUID


from amprenta_rag.config import get_config
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.ingestion.transcriptomics.embedding import embed_transcriptomics_dataset
from amprenta_rag.ingestion.transcriptomics.file_parsing import extract_gene_set_from_file
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_transcriptomics_dataset_page(
    file_path: str,
    gene_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal transcriptomics dataset.

    Args:
        file_path: Path to the transcriptomics file
        gene_count: Number of unique normalized genes
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Transcriptomics",
        entity_count=gene_count,
        raw_rows=raw_rows,
        entity_name="genes",
    )


def ingest_transcriptomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal transcriptomics dataset file.

    Args:
        file_path: Path to the transcriptomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Starting ingestion of transcriptomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Transcriptomics file not found: {file_path}")

    # Extract genes from file
    gene_set, df = extract_gene_set_from_file(file_path)

    if not gene_set:
        raise ValueError(f"No genes extracted from file: {file_path}")

    # Find gene column name for embedding
    gene_column = None
    candidate_names = ["gene", "Gene", "gene_name", "Gene Symbol", "GeneSymbol", "symbol", "Symbol"]
    for col in df.columns:
        if col in candidate_names or col.lower() in [c.lower() for c in candidate_names]:
            gene_column = col
            break
    if not gene_column:
        # Fallback: use first column
        gene_column = df.columns[0] if len(df.columns) > 0 else "gene"

    # TIER 3: Create dataset in Postgres FIRST (primary)
    cfg = get_config()
    postgres_dataset = None
    page_id: Optional[str] = None

    if cfg.pipeline.use_postgres_as_sot:
        try:
            from amprenta_rag.ingestion.postgres_integration import (
                create_or_update_dataset_in_postgres,
            )

            dataset_name = f"Internal Transcriptomics — {Path(file_path).stem}"
            postgres_dataset = create_or_update_dataset_in_postgres(
                name=dataset_name,
                omics_type=OmicsType.TRANSCRIPTOMICS,
                file_paths=[file_path],
                description=f"Internal transcriptomics dataset ingested from {file_path}. Contains {len(df)} raw entries, {len(gene_set)} normalized genes.",
                notion_page_id=notion_page_id,
                program_ids=None,
                experiment_ids=None,
            )
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Created/updated dataset in Postgres: %s (ID: %s)",
                postgres_dataset.name,
                postgres_dataset.id,
            )
            page_id = notion_page_id or str(postgres_dataset.id)
        except Exception as e:
            logger.error(
                "[INGEST][TRANSCRIPTOMICS] Postgres creation failed: %r",
                e,
            )
            if cfg.pipeline.use_postgres_as_sot:
                raise RuntimeError(f"Postgres creation failed (required): {e}") from e

    # Attach file (noted in summary for now) - only if Notion page exists
    if page_id and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        attach_file_to_page(page_id, file_path, "Transcriptomics")

    # Link to Programs and Experiments (Notion only if enabled)
    if (program_ids or experiment_ids) and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            link_to_programs_and_experiments(
                dataset_page_id=page_id,
                omics_type="Transcriptomics",
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )

    # Link to Programs/Experiments in Postgres (if Postgres dataset exists)
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset and (program_ids or experiment_ids):
        try:
            from amprenta_rag.ingestion.postgres_program_experiment_linking import (
                link_dataset_to_programs_and_experiments_in_postgres,
            )

            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Linking dataset %s to programs/experiments in Postgres",
                postgres_dataset.id,
            )

            link_results = link_dataset_to_programs_and_experiments_in_postgres(
                dataset_id=cast(UUID, postgres_dataset.id),
                notion_program_ids=program_ids,
                notion_experiment_ids=experiment_ids,
            )

            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Linked dataset to %d programs, %d experiments in Postgres",
                link_results["programs_linked"],
                link_results["experiments_linked"],
            )
        except Exception as e:
            logger.warning(
                "[INGEST][TRANSCRIPTOMICS] Postgres program/experiment linking skipped (error): %r",
                e,
            )

    # Link features to Postgres (if Postgres dataset exists)
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset and cfg.pipeline.enable_feature_linking:
        try:
            from amprenta_rag.ingestion.features.postgres_linking import (
                batch_link_features_to_dataset_in_postgres,
            )
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Batch linking %d genes to Postgres (dataset: %s)",
                len(gene_set),
                postgres_dataset.id,
            )

            feature_tuples = [(name, "gene") for name in gene_set]

            feature_results = batch_link_features_to_dataset_in_postgres(
                features=feature_tuples,
                dataset_id=cast(UUID, postgres_dataset.id),
                max_workers=cfg.pipeline.feature_linking_max_workers,
            )

            linked_count = sum(1 for v in feature_results.values() if v)
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Batch linked %d/%d genes to Postgres",
                linked_count,
                len(gene_set),
            )
        except Exception as e:
            logger.warning(
                "[INGEST][TRANSCRIPTOMICS] Postgres feature linking skipped (error): %r",
                e,
            )

    # Link genes to Gene Features DB in Notion (batch-optimized)
    # Only if Notion sync enabled (for backward compatibility)
    if cfg.pipeline.enable_feature_linking and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            try:
                from amprenta_rag.ingestion.features.batch_linking import batch_link_features

                logger.info(
                    "[INGEST][TRANSCRIPTOMICS] Batch linking %d genes to Notion Gene Features DB (max_workers=%d)",
                    len(gene_set),
                    cfg.pipeline.feature_linking_max_workers,
                )

                features = [("gene", gene) for gene in gene_set]
                feature_pages = batch_link_features(
                    features=features,
                    dataset_page_id=page_id,
                    max_workers=cfg.pipeline.feature_linking_max_workers,
                    enable_linking=True,
                )

                linked_count = sum(1 for v in feature_pages.values() if v)
                logger.info(
                    "[INGEST][TRANSCRIPTOMICS] Batch linked %d/%d genes to Notion Gene Features DB",
                    linked_count,
                    len(gene_set),
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][TRANSCRIPTOMICS] Notion feature linking skipped (error): %r",
                    e,
                )
        else:
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Notion feature linking skipped (no Notion page ID)"
            )
    elif not cfg.pipeline.enable_feature_linking:
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Feature linking disabled (ENABLE_FEATURE_LINKING=false)"
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem

    # TIER 3: Use Postgres-aware embedding if Postgres dataset exists
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset:
        try:
            from amprenta_rag.ingestion.postgres_integration import (
                embed_dataset_with_postgres_metadata,
            )
            embed_dataset_with_postgres_metadata(
                dataset_id=cast(UUID, postgres_dataset.id),
                dataset_name=dataset_name,
                species_or_features=list(gene_set),
                omics_type=OmicsType.TRANSCRIPTOMICS,
                signature_matches=None,  # TODO: Add signature matching for transcriptomics
                notion_page_id=page_id,
            )
            logger.info(
                "[INGEST][TRANSCRIPTOMICS] Embedded dataset to Pinecone using Postgres metadata"
            )
        except Exception as e:
            logger.warning(
                "[INGEST][TRANSCRIPTOMICS] Postgres embedding failed, falling back to Notion: %r",
                e,
            )
            # Fallback to Notion-based embedding
            embed_transcriptomics_dataset(
                page_id=page_id or "",
                dataset_name=dataset_name,
                genes=gene_set,
                df=df,
                gene_column=gene_column,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
    else:
        # Use Notion-based embedding (default)
        embed_transcriptomics_dataset(
            page_id=page_id or "",
            dataset_name=dataset_name,
            genes=gene_set,
            df=df,
            gene_column=gene_column,
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    logger.info(
        "[INGEST][TRANSCRIPTOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id or ""

