"""
Main metabolomics ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
metabolomics dataset ingestion pipeline: file parsing, Notion page creation,
feature linking, and RAG embedding.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import fetch_dataset_page
from amprenta_rag.ingestion.metabolomics.embedding import embed_metabolomics_dataset
from amprenta_rag.ingestion.metabolomics.file_parsing import extract_metabolite_set_from_file
from amprenta_rag.ingestion.omics_ingestion_utils import (
    attach_file_to_page,
    create_omics_dataset_page,
    link_to_programs_and_experiments,
)
from amprenta_rag.logging_utils import get_logger

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

    Args:
        file_path: Path to the metabolomics file (CSV/TSV)
        notion_page_id: Existing Experimental Data Assets page ID to attach to
        create_page: If True and notion_page_id is None, create a new page
        program_ids: Optional list of Program page IDs to link
        experiment_ids: Optional list of Experiment page IDs to link

    Returns:
        The Notion page ID of the Experimental Data Asset
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

    # Determine page ID
    if notion_page_id:
        # Use existing page
        page_id = notion_page_id
        logger.info(
            "[INGEST][METABOLOMICS] Using existing dataset page %s",
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
                f"Added metabolomics file {file_basename} with {len(metabolite_set)} unique metabolites."
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
                "[INGEST][METABOLOMICS] Could not update summary for page %s: %r",
                page_id,
                e,
            )
    elif create_page:
        # Create new page
        page_id = create_metabolomics_dataset_page(
            file_path=file_path,
            metabolite_count=len(metabolite_set),
            raw_rows=raw_rows,
        )
    else:
        raise ValueError(
            "Either notion_page_id must be provided or create_page must be True"
        )

    # Attach file (noted in summary for now)
    attach_file_to_page(page_id, file_path, "Metabolomics")

    # Link to Programs and Experiments
    if program_ids or experiment_ids:
        link_to_programs_and_experiments(
            dataset_page_id=page_id,
            omics_type="Metabolomics",
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

    # Link metabolites to Metabolite Features DB
    try:
        from amprenta_rag.ingestion.feature_extraction import link_feature

        logger.info(
            "[INGEST][METABOLOMICS] Linking %d metabolites to Metabolite Features DB",
            len(metabolite_set),
        )
        linked_count = 0
        for metabolite in metabolite_set:
            try:
                link_feature("metabolite", metabolite, page_id)
                linked_count += 1
            except Exception as e:
                logger.warning(
                    "[INGEST][METABOLOMICS] Error linking metabolite '%s': %r",
                    metabolite,
                    e,
                )
        logger.info(
            "[INGEST][METABOLOMICS] Linked %d/%d metabolites to Metabolite Features DB",
            linked_count,
            len(metabolite_set),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][METABOLOMICS] Feature linking skipped (error): %r",
            e,
        )

    # Embed into Pinecone
    dataset_name = Path(file_path).stem
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

