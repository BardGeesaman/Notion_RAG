"""
Main lipidomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
lipidomics dataset ingestion pipeline: file parsing, Postgres storage,
feature linking, signature scoring, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.config import get_config
from amprenta_rag.ingestion.lipidomics.embedding import embed_lipidomics_dataset
from amprenta_rag.ingestion.lipidomics.file_parsing import extract_species_from_file
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset,
    update_dataset_with_signature_matches,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_lipidomics_dataset_page(
    file_path: str,
    species_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal lipidomics dataset.

    Args:
        file_path: Path to the lipidomics file
        species_count: Number of unique normalized species
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Lipidomics",
        entity_count=species_count,
        raw_rows=raw_rows,
        entity_name="species",
    )


def ingest_lipidomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal lipidomics dataset file.

    Args:
        file_path: Path to the lipidomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][LIPIDOMICS] Starting ingestion of lipidomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Lipidomics file not found: {file_path}")

    # Extract species from file
    species_set, raw_rows = extract_species_from_file(file_path)

    if not species_set:
        raise ValueError(f"No species extracted from file: {file_path}")

    # TIER 3: Create dataset in Postgres FIRST (primary)
    cfg = get_config()
    postgres_dataset = None
    page_id = None
    
    if cfg.pipeline.use_postgres_as_sot:
        try:
            from amprenta_rag.ingestion.postgres_integration import (
                create_or_update_dataset_in_postgres,
            )
            
            dataset_name = f"Internal Lipidomics — {Path(file_path).stem}"
            postgres_dataset = create_or_update_dataset_in_postgres(
                name=dataset_name,
                omics_type=OmicsType.LIPIDOMICS,
                file_paths=[file_path],
                description=f"Internal lipidomics dataset ingested from {file_path}. Contains {raw_rows} raw entries, {len(species_set)} normalized species.",
                notion_page_id=notion_page_id,  # Link to existing Notion if provided
                program_ids=None,
                experiment_ids=None,
            )
            logger.info(
                "[INGEST][LIPIDOMICS] Created/updated dataset in Postgres: %s (ID: %s)",
                postgres_dataset.name,
                postgres_dataset.id,
            )
            # Use Postgres ID as primary identifier (fallback to Notion ID if provided)
            page_id = notion_page_id or str(postgres_dataset.id)
        except Exception as e:
            logger.error(
                "[INGEST][LIPIDOMICS] Postgres creation failed: %r",
                e,
            )
            # If Postgres is required, fail; otherwise fallback to Notion
            if cfg.pipeline.use_postgres_as_sot:
                raise RuntimeError(f"Postgres creation failed (required): {e}") from e
    
    # Create/update Notion page (OPTIONAL - only if enabled)
    if cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write:
        try:
            if notion_page_id:
                # Use existing page
                page_id = notion_page_id
                logger.info(
                    "[INGEST][LIPIDOMICS] Using existing Notion dataset page %s",
                    page_id,
                )
                # Update summary
                try:
                    page = fetch_dataset_page(page_id)
                    page_props = page.get("properties", {}) or {}
                    summary_prop = page_props.get("Summary", {}) or {}
                    existing_summary = "".join(
                        rt.get("plain_text", "")
                        for rt in summary_prop.get("rich_text", [])
                    )
                    new_note = (
                        f"\n\n---\n\n"
                        f"Lipidomics file ingested: {file_path}\n"
                        f"Contains {raw_rows} raw entries, {len(species_set)} normalized species."
                    )
                    combined_summary = existing_summary + new_note

                    # Notion update skipped - Notion support removed
                    logger.debug("[INGEST][LIPIDOMICS] Notion summary update skipped (Notion removed)")
                except Exception as e:
                    logger.warning(
                        "[INGEST][LIPIDOMICS] Could not update Notion summary for page %s: %r",
                        page_id,
                        e,
                    )
            elif create_page:
                # Create new Notion page
                page_id = create_lipidomics_dataset_page(
                    file_path=file_path,
                    species_count=len(species_set),
                    raw_rows=raw_rows,
                )
                # Update Postgres with Notion page ID if Postgres dataset exists
                if postgres_dataset:
                    from amprenta_rag.database.session import db_session
                    with db_session() as db:
                        postgres_dataset.notion_page_id = page_id
                        db.commit()
                        logger.info(
                            "[INGEST][LIPIDOMICS] Linked Postgres dataset %s to Notion page %s",
                            postgres_dataset.id,
                            page_id,
                        )
            elif not cfg.pipeline.use_postgres_as_sot:
                # Notion-only mode: require page creation
                raise ValueError(
                    "Either notion_page_id must be provided or create_page must be True"
                )
        except Exception as e:
            logger.warning(
                "[INGEST][LIPIDOMICS] Notion sync skipped (error): %r",
                e,
            )
            # Non-blocking - continue with Postgres-only mode
    elif not cfg.pipeline.use_postgres_as_sot:
        # Notion-only mode (legacy)
        if notion_page_id:
            page_id = notion_page_id
        elif create_page:
            page_id = create_lipidomics_dataset_page(
                file_path=file_path,
                species_count=len(species_set),
                raw_rows=raw_rows,
            )
        else:
            raise ValueError(
                "Either notion_page_id must be provided or create_page must be True"
            )
    
    # Attach file (noted in summary for now) - only if Notion page exists
    if page_id and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        attach_file_to_page(page_id, file_path, "Lipidomics")

    # Link to Programs and Experiments (Notion only if enabled)
    if (program_ids or experiment_ids) and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            link_to_programs_and_experiments(
                dataset_page_id=page_id,
                omics_type="Lipidomics",
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
                "[INGEST][LIPIDOMICS] Linking dataset %s to programs/experiments in Postgres",
                postgres_dataset.id,
            )
            
            results = link_dataset_to_programs_and_experiments_in_postgres(
                dataset_id=postgres_dataset.id,
                notion_program_ids=program_ids,
                notion_experiment_ids=experiment_ids,
            )
            
            logger.info(
                "[INGEST][LIPIDOMICS] Linked dataset to %d programs, %d experiments in Postgres",
                results["programs_linked"],
                results["experiments_linked"],
            )
        except Exception as e:
            logger.warning(
                "[INGEST][LIPIDOMICS] Postgres program/experiment linking skipped (error): %r",
                e,
            )

    # Link features to Postgres (if Postgres dataset exists)
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset and cfg.pipeline.enable_feature_linking:
        try:
            from amprenta_rag.ingestion.features.postgres_linking import (
                batch_link_features_to_dataset_in_postgres,
            )
            from amprenta_rag.models.domain import FeatureType
            
            logger.info(
                "[INGEST][LIPIDOMICS] Batch linking %d lipid species to Postgres (dataset: %s)",
                len(species_set),
                postgres_dataset.id,
            )
            
            feature_tuples = [(name, FeatureType.LIPID) for name in species_set]
            
            results = batch_link_features_to_dataset_in_postgres(
                features=feature_tuples,
                dataset_id=postgres_dataset.id,
                max_workers=cfg.pipeline.feature_linking_max_workers,
            )
            
            linked_count = sum(1 for v in results.values() if v)
            logger.info(
                "[INGEST][LIPIDOMICS] Batch linked %d/%d lipid species to Postgres",
                linked_count,
                len(species_set),
            )
        except Exception as e:
            logger.warning(
                "[INGEST][LIPIDOMICS] Postgres feature linking skipped (error): %r",
                e,
            )

    # Link lipid species to Lipid Species DB in Notion (batch-optimized)
    # Only if Notion sync enabled (for backward compatibility)
    if cfg.pipeline.enable_feature_linking and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            try:
                from amprenta_rag.ingestion.features.batch_linking import batch_link_features

                logger.info(
                    "[INGEST][LIPIDOMICS] Batch linking %d lipid species to Notion Lipid Species DB (max_workers=%d)",
                    len(species_set),
                    cfg.pipeline.feature_linking_max_workers,
                )
                
                features = [("lipid", lipid) for lipid in species_set]
                feature_pages = batch_link_features(
                    features=features,
                    dataset_page_id=page_id,
                    max_workers=cfg.pipeline.feature_linking_max_workers,
                    enable_linking=True,
                )
                
                linked_count = sum(1 for v in feature_pages.values() if v)
                logger.info(
                    "[INGEST][LIPIDOMICS] Batch linked %d/%d lipid species to Notion Lipid Species DB",
                    linked_count,
                    len(species_set),
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][LIPIDOMICS] Notion feature linking skipped (error): %r",
                    e,
                )
        else:
            logger.info(
                "[INGEST][LIPIDOMICS] Notion feature linking skipped (no Notion page ID)"
            )
    elif not cfg.pipeline.enable_feature_linking:
        logger.info(
            "[INGEST][LIPIDOMICS] Feature linking disabled (ENABLE_FEATURE_LINKING=false)"
        )

    # Score against signatures
    signature_matches = []
    if cfg.pipeline.enable_signature_scoring:
        # Use Postgres dataset ID if available, otherwise Notion page ID
        dataset_identifier = str(postgres_dataset.id) if postgres_dataset else page_id
        logger.info(
            "[INGEST][LIPIDOMICS] Scoring dataset %s against signatures (%d species)",
            dataset_identifier,
            len(species_set),
        )

        matches = find_matching_signatures_for_dataset(
            dataset_species=species_set,  # Legacy support
            overlap_threshold=cfg.pipeline.signature_overlap_threshold,
            dataset_page_id=page_id if page_id else dataset_identifier,  # Multi-omics support
            omics_type="Lipidomics",
        )

        if matches:
            logger.info(
                "[INGEST][LIPIDOMICS] Found %d matching signature(s) for dataset %s",
                len(matches),
                dataset_identifier,
            )
            # Update Notion if enabled
            if page_id and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
                update_dataset_with_signature_matches(
                    dataset_page_id=page_id,
                    matches=matches,
                )
            # TODO: Update Postgres dataset with signature matches
            signature_matches = matches
        else:
            logger.info(
                "[INGEST][LIPIDOMICS] No signature matches found for dataset %s",
                dataset_identifier,
            )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    
    # TIER 3: Use Postgres-aware embedding if Postgres dataset exists
    cfg = get_config()
    if cfg.pipeline.use_postgres_as_sot and postgres_dataset:
        try:
            from amprenta_rag.ingestion.postgres_integration import (
                embed_dataset_with_postgres_metadata,
            )
            embed_dataset_with_postgres_metadata(
                dataset_id=postgres_dataset.id,
                dataset_name=dataset_name,
                species_or_features=list(species_set),
                omics_type=OmicsType.LIPIDOMICS,
                signature_matches=signature_matches,
                notion_page_id=page_id,
            )
            logger.info(
                "[INGEST][LIPIDOMICS] Embedded dataset to Pinecone using Postgres metadata"
            )
        except Exception as e:
            logger.warning(
                "[INGEST][LIPIDOMICS] Postgres embedding failed, falling back to Notion: %r",
                e,
            )
            # Fallback to Notion-based embedding
            embed_lipidomics_dataset(
                page_id=page_id,
                dataset_name=dataset_name,
                species=species_set,
                signature_matches=signature_matches,
            )
    else:
        # Use Notion-based embedding (default)
        embed_lipidomics_dataset(
            page_id=page_id,
            dataset_name=dataset_name,
            species=species_set,
            signature_matches=signature_matches,
        )

    # Auto-ingest linked Experiment pages (only if Notion sync enabled)
    if experiment_ids and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        from amprenta_rag.ingestion.experiments_ingestion import ingest_experiment

        for exp_id in experiment_ids:
            try:
                logger.info(
                    "[INGEST][EXPERIMENT] Auto-ingesting Experiment %s after lipidomics dataset %s created",
                    exp_id,
                    page_id,
                )
                ingest_experiment(
                    exp_page_id=exp_id,
                    parent_type="Experiment",
                    force=False,
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][EXPERIMENT] Error auto-ingesting Experiment %s: %r",
                    exp_id,
                    e,
                )
                # Non-blocking - continue with other experiments

    logger.info(
        "[INGEST][LIPIDOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

