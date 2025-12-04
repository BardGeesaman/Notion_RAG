"""
Main lipidomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
lipidomics dataset ingestion pipeline: file parsing, Notion page creation,
feature linking, signature scoring, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page
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

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][LIPIDOMICS] Using existing dataset page %s",
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
                "[INGEST][LIPIDOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_lipidomics_dataset_page(
            file_path=file_path,
            species_count=len(species_set),
            raw_rows=raw_rows,
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Lipidomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Lipidomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Link lipid species to Lipid Species DB
    try:
        from amprenta_rag.ingestion.feature_extraction import link_feature

        logger.info(
            "[INGEST][LIPIDOMICS] Linking %d lipid species to Lipid Species DB",
            len(species_set),
        )
        linked_count = 0
        for lipid in species_set:
            try:
                link_feature("lipid", lipid, page_id)
                linked_count += 1
            except Exception as e:
                logger.warning(
                    "[INGEST][LIPIDOMICS] Error linking lipid species '%s': %r",
                    lipid,
                    e,
                )
        logger.info(
            "[INGEST][LIPIDOMICS] Linked %d/%d lipid species to Lipid Species DB",
            linked_count,
            len(species_set),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][LIPIDOMICS] Feature linking skipped (error): %r",
            e,
        )

    # Score against signatures
    cfg = get_config()
    signature_matches = []
    if cfg.pipeline.enable_signature_scoring:
        logger.info(
            "[INGEST][LIPIDOMICS] Scoring dataset %s against signatures (%d species)",
            page_id,
            len(species_set),
        )

        matches = find_matching_signatures_for_dataset(
            dataset_species=species_set,  # Legacy support
            overlap_threshold=cfg.pipeline.signature_overlap_threshold,
            dataset_page_id=page_id,  # Multi-omics support
            omics_type="Lipidomics",
        )

        if matches:
            logger.info(
                "[INGEST][LIPIDOMICS] Found %d matching signature(s) for dataset %s",
                len(matches),
                page_id,
            )
            update_dataset_with_signature_matches(
                dataset_page_id=page_id,
                matches=matches,
            )
            signature_matches = matches
        else:
            logger.info(
                "[INGEST][LIPIDOMICS] No signature matches found for dataset %s",
                page_id,
            )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
    embed_lipidomics_dataset(
        page_id=page_id,
        dataset_name=dataset_name,
        species=species_set,
        signature_matches=signature_matches,
    )

    # Auto-ingest linked Experiment pages
    if experiment_ids:
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

