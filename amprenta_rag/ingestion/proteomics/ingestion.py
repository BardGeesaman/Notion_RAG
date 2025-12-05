"""
Main proteomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
proteomics dataset ingestion pipeline: file parsing, Notion page creation,
feature linking, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.ingestion.proteomics.embedding import embed_proteomics_dataset
from amprenta_rag.ingestion.proteomics.file_parsing import extract_protein_set_from_file
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)


def create_proteomics_dataset_page(
    file_path: str,
    protein_count: int,
    raw_rows: int,
) -> str:
    """
    Create a new Experimental Data Asset page for an internal proteomics dataset.

    Args:
        file_path: Path to the proteomics file
        protein_count: Number of unique normalized proteins
        raw_rows: Number of raw rows in the file

    Returns:
        Created Notion page ID
    """
    return create_omics_dataset_page(
        file_path=file_path,
        omics_type="Proteomics",
        entity_count=protein_count,
        raw_rows=raw_rows,
        entity_name="proteins",
    )


def ingest_proteomics_file(
    file_path: str,
    notion_page_id: Optional[str] = None,
    create_page: bool = False,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> str:
    """
    Ingest an internal proteomics dataset file.

    Args:
        file_path: Path to the proteomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
    """
    logger.info(
        "[INGEST][PROTEOMICS] Starting ingestion of proteomics file: %s",
        file_path,
    )

    # Validate file exists
    file_path_obj = Path(file_path)
    if not file_path_obj.exists():
        raise FileNotFoundError(f"Proteomics file not found: {file_path}")

    # Extract proteins from file
    protein_set, raw_rows = extract_protein_set_from_file(file_path)

    if not protein_set:
        raise ValueError(f"No proteins extracted from file: {file_path}")

    # TIER 3: Create dataset in Postgres FIRST (primary)
    cfg = get_config()
    postgres_dataset = None
    page_id = None
    
    if cfg.pipeline.use_postgres_as_sot:
        try:
            from amprenta_rag.ingestion.postgres_integration import (
                create_or_update_dataset_in_postgres,
            )
            
            dataset_name = f"Internal Proteomics — {Path(file_path).stem}"
            postgres_dataset = create_or_update_dataset_in_postgres(
                name=dataset_name,
                omics_type=OmicsType.PROTEOMICS,
                file_paths=[file_path],
                description=f"Internal proteomics dataset ingested from {file_path}. Contains {raw_rows} raw entries, {len(protein_set)} normalized proteins.",
                notion_page_id=notion_page_id,
                program_ids=None,  # TODO: Convert Notion IDs to Postgres UUIDs
                experiment_ids=None,
            )
            logger.info(
                "[INGEST][PROTEOMICS] Created/updated dataset in Postgres: %s (ID: %s)",
                postgres_dataset.name,
                postgres_dataset.id,
            )
            page_id = notion_page_id or str(postgres_dataset.id)
        except Exception as e:
            logger.error(
                "[INGEST][PROTEOMICS] Postgres creation failed: %r",
                e,
            )
            if cfg.pipeline.use_postgres_as_sot:
                raise RuntimeError(f"Postgres creation failed (required): {e}") from e

    # Create/update Notion page (OPTIONAL - only if enabled)
    if cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write:
        try:
            if notion_page_id:
                # Use existing page
                page_id = notion_page_id
                logger.info(
                    "[INGEST][PROTEOMICS] Using existing Notion dataset page %s",
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
                    file_basename = Path(file_path).name
                    new_note = (
                        f"\n\n---\n\n"
                        f"Added proteomics file {file_basename} with {len(protein_set)} unique proteins."
                    )
                    combined_summary = existing_summary + new_note

                    # Update page
                    url = f"{get_config().notion.base_url}/pages/{page_id}"
                    resp = requests.patch(
                        url,
                        headers=notion_headers(),
                        json={"properties": {"Summary": {"rich_text": [{"text": {"content": combined_summary}}]}}},
                        timeout=30,
                    )
                    resp.raise_for_status()
                except Exception as e:
                    logger.warning(
                        "[INGEST][PROTEOMICS] Could not update Notion summary for page %s: %r",
                        page_id,
                        e,
                    )
            elif create_page:
                # Create new Notion page
                page_id = create_proteomics_dataset_page(
                    file_path=file_path,
                    protein_count=len(protein_set),
                    raw_rows=raw_rows,
                )
                # Update Postgres with Notion page ID if Postgres dataset exists
                if postgres_dataset:
                    from amprenta_rag.database.base import get_db
                    db = next(get_db())
                    postgres_dataset.notion_page_id = page_id
                    db.commit()
                    logger.info(
                        "[INGEST][PROTEOMICS] Linked Postgres dataset %s to Notion page %s",
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
                "[INGEST][PROTEOMICS] Notion sync skipped (error): %r",
                e,
            )
    elif not cfg.pipeline.use_postgres_as_sot:
        # Notion-only mode (legacy)
        if notion_page_id:
            page_id = notion_page_id
        elif create_page:
            page_id = create_proteomics_dataset_page(
                file_path=file_path,
                protein_count=len(protein_set),
                raw_rows=raw_rows,
            )
        else:
            raise ValueError(
                "Either notion_page_id must be provided or create_page must be True"
            )
    
    # Attach file (noted in summary for now) - only if Notion page exists
    if page_id and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        attach_file_to_page(page_id, file_path, "Proteomics")

    # Link to Programs and Experiments (Notion only if enabled)
    if (program_ids or experiment_ids) and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            link_to_programs_and_experiments(
                dataset_page_id=page_id,
                omics_type="Proteomics",
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
                "[INGEST][PROTEOMICS] Linking dataset %s to programs/experiments in Postgres",
                postgres_dataset.id,
            )
            
            results = link_dataset_to_programs_and_experiments_in_postgres(
                dataset_id=postgres_dataset.id,
                notion_program_ids=program_ids,
                notion_experiment_ids=experiment_ids,
            )
            
            logger.info(
                "[INGEST][PROTEOMICS] Linked dataset to %d programs, %d experiments in Postgres",
                results["programs_linked"],
                results["experiments_linked"],
            )
        except Exception as e:
            logger.warning(
                "[INGEST][PROTEOMICS] Postgres program/experiment linking skipped (error): %r",
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
                "[INGEST][PROTEOMICS] Batch linking %d proteins to Postgres (dataset: %s)",
                len(protein_set),
                postgres_dataset.id,
            )
            
            feature_tuples = [(name, FeatureType.PROTEIN) for name in protein_set]
            
            results = batch_link_features_to_dataset_in_postgres(
                features=feature_tuples,
                dataset_id=postgres_dataset.id,
                max_workers=cfg.pipeline.feature_linking_max_workers,
            )
            
            linked_count = sum(1 for v in results.values() if v)
            logger.info(
                "[INGEST][PROTEOMICS] Batch linked %d/%d proteins to Postgres",
                linked_count,
                len(protein_set),
            )
        except Exception as e:
            logger.warning(
                "[INGEST][PROTEOMICS] Postgres feature linking skipped (error): %r",
                e,
            )

    # Link proteins to Protein Features DB in Notion (batch-optimized)
    # Only if Notion sync enabled (for backward compatibility)
    if cfg.pipeline.enable_feature_linking and (cfg.pipeline.enable_notion_sync or cfg.pipeline.enable_dual_write or not cfg.pipeline.use_postgres_as_sot):
        if page_id:
            try:
                from amprenta_rag.ingestion.features.batch_linking import batch_link_features

                logger.info(
                    "[INGEST][PROTEOMICS] Batch linking %d proteins to Notion Protein Features DB (max_workers=%d)",
                    len(protein_set),
                    cfg.pipeline.feature_linking_max_workers,
                )
                
                features = [("protein", protein) for protein in protein_set]
                feature_pages = batch_link_features(
                    features=features,
                    dataset_page_id=page_id,
                    max_workers=cfg.pipeline.feature_linking_max_workers,
                    enable_linking=True,
                )
                
                linked_count = sum(1 for v in feature_pages.values() if v)
                logger.info(
                    "[INGEST][PROTEOMICS] Batch linked %d/%d proteins to Notion Protein Features DB",
                    linked_count,
                    len(protein_set),
                )
            except Exception as e:
                logger.warning(
                    "[INGEST][PROTEOMICS] Notion feature linking skipped (error): %r",
                    e,
                )
        else:
            logger.info(
                "[INGEST][PROTEOMICS] Notion feature linking skipped (no Notion page ID)"
            )
    elif not cfg.pipeline.enable_feature_linking:
        logger.info(
            "[INGEST][PROTEOMICS] Feature linking disabled (ENABLE_FEATURE_LINKING=false)"
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
                dataset_id=postgres_dataset.id,
                dataset_name=dataset_name,
                species_or_features=list(protein_set),
                omics_type=OmicsType.PROTEOMICS,
                signature_matches=None,  # TODO: Add signature matching for proteomics
                notion_page_id=page_id,
            )
            logger.info(
                "[INGEST][PROTEOMICS] Embedded dataset to Pinecone using Postgres metadata"
            )
        except Exception as e:
            logger.warning(
                "[INGEST][PROTEOMICS] Postgres embedding failed, falling back to Notion: %r",
                e,
            )
            # Fallback to Notion-based embedding
            embed_proteomics_dataset(
                page_id=page_id,
                dataset_name=dataset_name,
                proteins=protein_set,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
    else:
        # Use Notion-based embedding (default)
        embed_proteomics_dataset(
            page_id=page_id,
            dataset_name=dataset_name,
            proteins=protein_set,
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    logger.info(
        "[INGEST][PROTEOMICS] Completed ingestion of file %s -> dataset page %s",
        file_path,
        page_id,
    )

    return page_id

