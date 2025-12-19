"""
Main metabolomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
metabolomics dataset ingestion pipeline: file parsing, Postgres storage,
feature linking, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional


from amprenta_rag.config import get_config
from amprenta_rag.ingestion.metabolomics.embedding import embed_metabolomics_dataset
from amprenta_rag.ingestion.metabolomics.file_parsing import extract_metabolite_set_from_file
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.ingestion.postgres_integration import (
    create_or_update_dataset_in_postgres,
    embed_dataset_with_postgres_metadata,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_metabolomics_dataset_page(
    file_path: str,
    metabolite_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal metabolomics dataset.

    Args:
        file_path: Path to the metabolomics file
        metabolite_count: Number of unique normalized metabolites
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Metabolomics",
        entity_count=metabolite_count,
        raw_rows=raw_rows,
        entity_name="metabolites",
    )


def ingest_metabolomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal polar metabolomics dataset file.

    This function orchestrates the complete metabolomics ingestion pipeline:
    1. Validates and parses the input file
    2. Extracts and normalizes metabolite names
    3. Creates/updates Postgres dataset (if USE_POSTGRES_AS_SOT=true)
    4. Optionally creates/updates Notion page (if ENABLE_NOTION_SYNC=true)
    5. Links features to Notion feature databases (if enabled)
    6. Embeds dataset into Pinecone for RAG

    Architecture:
    - Postgres is the primary database (if USE_POSTGRES_AS_SOT=true)
    - Notion is optional for documentation (if ENABLE_NOTION_SYNC=true)
    - Feature linking uses Notion (future: will use Postgres)

    Args:
        file_path: Path to the metabolomics file (CSV/TSV). Must exist and be readable.
        notion_page_id: Optional existing Notion Experimental Data Assets page ID to link to.
            If provided, the dataset will be linked to this page.
        create_page: If True and notion_page_id is None, create a new Notion page.
            Requires ENABLE_NOTION_SYNC=true. Default: False.
        program_ids: Optional list of Notion Program page IDs to link to the dataset.
            Requires Notion sync enabled. Default: None.
        experiment_ids: Optional list of Notion Experiment page IDs to link to the dataset.
            Requires Notion sync enabled. Default: None.

    Returns:
        The dataset identifier (Notion page ID if Notion sync enabled, Postgres UUID otherwise).
        For Postgres-only mode, returns the Postgres dataset UUID as a string.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If no metabolites can be extracted from the file.
        RuntimeError: If Postgres creation fails and USE_POSTGRES_AS_SOT=true.

    Example:
        >>> # Postgres-only ingestion
        >>> dataset_id = ingest_metabolomics_file(
        ...     file_path="data/metabolomics_sample.csv",
        ...     create_page=False,
        ... )

        >>> # With Notion sync (requires ENABLE_NOTION_SYNC=true)
        >>> dataset_id = ingest_metabolomics_file(
        ...     file_path="data/metabolomics_sample.csv",
        ...     create_page=True,
        ...     program_ids=["program-page-id"],
        ... )
    """
    logger.info(
        "[INGEST][METABOLOMICS] Starting ingestion of metabolomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Metabolomics file not found: {file_path}")

    # Extract metabolites from file
    metabolite_set, raw_rows = extract_metabolite_set_from_file(file_path)

    if not metabolite_set:
        raise ValueError(f"No metabolites extracted from file: {file_path}")

    # TIER 3: Create dataset in Postgres FIRST (primary)
    cfg = get_config()
    postgres_dataset = None
    page_id = None

    if cfg.pipeline.use_postgres_as_sot:
        try:
            dataset_name = f"Internal Metabolomics — {Path(file_path).stem}"
            postgres_dataset = create_or_update_dataset_in_postgres(
                name=dataset_name,
                omics_type=OmicsType.METABOLOMICS,
                file_paths=[file_path],
                description=f"Internal metabolomics dataset ingested from {file_path}. Contains {raw_rows} raw entries, {len(metabolite_set)} normalized metabolites.",
                notion_page_id=notion_page_id,
                program_ids=None,
                experiment_ids=None,
            )
            logger.info(
                "[INGEST][METABOLOMICS] Created/updated dataset in Postgres: %s (ID: %s)",
                postgres_dataset.name,
                postgres_dataset.id,
            )
            page_id = notion_page_id or str(postgres_dataset.id)
        except Exception as e:
            logger.error(
                "[INGEST][METABOLOMICS] Postgres creation failed: %r",
                e,
            )
            if cfg.pipeline.use_postgres_as_sot:
                raise RuntimeError(f"Postgres creation failed (required): {e}") from e

    # Attach file (noted in summary for now) - only if Notion page exists
    if page_id and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        attach_file_to_page(page_id, file_path, "Metabolomics")

    # Link to Programs and Experiments (Notion only if enabled)
    if (program_ids or experiment_ids) and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            link_to_programs_and_experiments(
                dataset_page_id=page_id,
                omics_type="Metabolomics",
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
                "[INGEST][METABOLOMICS] Linking dataset %s to programs/experiments in Postgres",
                postgres_dataset.id,
            )

            results = link_dataset_to_programs_and_experiments_in_postgres(
                dataset_id=postgres_dataset.id,
                notion_program_ids=program_ids,
                notion_experiment_ids=experiment_ids,
            )

            logger.info(
                "[INGEST][METABOLOMICS] Linked dataset to %d programs, %d experiments in Postgres",
                results["programs_linked"],
                results["experiments_linked"],
            )
        except Exception as e:
            logger.warning(
                "[INGEST][METABOLOMICS] Postgres program/experiment linking skipped (error): %r",
                e,
            )

    # Link features to Postgres (if Postgres dataset exists)
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset and cfg.pipeline.enable_feature_linking:
        try:
            from amprenta_rag.ingestion.features.postgres_linking import (
                batch_link_features_to_dataset_in_postgres,
            )

            logger.info(
                "[INGEST][METABOLOMICS] Batch linking %d metabolites to Postgres (dataset: %s)",
                len(metabolite_set),
                postgres_dataset.id,
            )

            # Convert to FeatureType tuples
            from amprenta_rag.models.domain import FeatureType
            feature_tuples = [(name, FeatureType.METABOLITE) for name in metabolite_set]

            results = batch_link_features_to_dataset_in_postgres(
                features=feature_tuples,
                dataset_id=postgres_dataset.id,
                max_workers=cfg.pipeline.feature_linking_max_workers,
            )

            linked_count = sum(1 for v in results.values() if v)
            logger.info(
                "[INGEST][METABOLOMICS] Batch linked %d/%d metabolites to Postgres",
                linked_count,
                len(metabolite_set),
            )
        except Exception as e:
            logger.warning(
                "[INGEST][METABOLOMICS] Postgres feature linking skipped (error): %r",
                e,
            )

    # Link metabolites to Metabolite Features DB in Notion (batch-optimized)
    # Only if Notion sync enabled (for backward compatibility)
    if cfg.pipeline.enable_feature_linking and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            try:
                from amprenta_rag.ingestion.features.batch_linking import batch_link_features

                logger.info(
                    "[INGEST][METABOLOMICS] Batch linking %d metabolites to Notion Metabolite Features DB (max_workers=%d)",
                    len(metabolite_set),
                    cfg.pipeline.feature_linking_max_workers,
                )

                features = [("metabolite", metabolite) for metabolite in metabolite_set]
                feature_pages = batch_link_features(
                    features=features,
                    dataset_page_id=page_id,
                    max_workers=cfg.pipeline.feature_linking_max_workers,
                    enable_linking=True,
                )

                linked_count = sum(1 for v in feature_pages.values() if v)
                logger.info(
                    "[INGEST][METABOLOMICS] Batch linked %d/%d metabolites to Notion Metabolite Features DB",
                    linked_count,
                    len(metabolite_set),
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][METABOLOMICS] Notion feature linking skipped (error): %r",
                    e,
                )
        else:
            logger.info(
                "[INGEST][METABOLOMICS] Notion feature linking skipped (no Notion page ID)"
            )
    elif not cfg.pipeline.enable_feature_linking:
        logger.info(
            "[INGEST][METABOLOMICS] Feature linking disabled (ENABLE_FEATURE_LINKING=false)"
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem

    # TIER 3: Use Postgres-aware embedding if Postgres dataset exists
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset:
        try:
            embed_dataset_with_postgres_metadata(
                dataset_id=postgres_dataset.id,
                dataset_name=dataset_name,
                species_or_features=list(metabolite_set),
                omics_type=OmicsType.METABOLOMICS,
                signature_matches=None,  # TODO: Add signature matching for metabolomics
                notion_page_id=page_id,
            )
            logger.info(
                "[INGEST][METABOLOMICS] Embedded dataset to Pinecone using Postgres metadata"
            )
        except Exception as e:
            logger.warning(
                "[INGEST][METABOLOMICS] Postgres embedding failed, falling back to Notion: %r",
                e,
            )
            # Fallback to Notion-based embedding
            embed_metabolomics_dataset(
                page_id=page_id,
                dataset_name=dataset_name,
                metabolites=metabolite_set,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
    else:
        # Use Notion-based embedding (default)
        embed_metabolomics_dataset(
            page_id=page_id,
            dataset_name=dataset_name,
            metabolites=metabolite_set,
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    logger.info(
        "[INGEST][METABOLOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

