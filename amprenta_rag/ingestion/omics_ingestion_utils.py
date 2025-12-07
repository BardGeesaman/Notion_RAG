"""
Shared utilities for omics dataset ingestion pipelines.

This module provides common functions used across lipidomics, metabolomics,
proteomics, and transcriptomics ingestion modules to reduce code duplication.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[OMICS-INGESTION-UTILS] notion_headers() deprecated - Notion support removed")
    return {}

__all__ = [
    "create_omics_dataset_page",
    "attach_file_to_page",
    "link_to_programs_and_experiments",
]


def create_omics_dataset_page(
    file_path: str,
    omics_type: str,
    entity_count: int,
    raw_rows: int,
    entity_name: str = "entities",
) -> str:
    """
    Create a new Experimental Data Asset page for an internal omics dataset.

    Args:
        file_path: Path to the omics file
        omics_type: Type of omics (e.g., "Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics")
        entity_count: Number of unique normalized entities (lipids, metabolites, proteins, genes)
        raw_rows: Number of raw rows in the file
        entity_name: Name of the entity type (e.g., "lipids", "metabolites", "proteins", "genes")

    Returns:
        Created Notion page ID
    """
    cfg = get_config().notion
    if not cfg.exp_data_db_id:
        raise RuntimeError("Experimental Data Assets DB ID not configured")

    file_name = Path(file_path).stem
    dataset_name = f"Internal {omics_type} — {file_name}"

    # Build properties
    props: Dict[str, Any] = {
        "Experiment Name": {"title": [{"text": {"content": dataset_name}}]},
        "Data Origin": {"select": {"name": "Internal – Amprenta"}},
        "Dataset Source Type": {"select": {"name": "Processed table"}},
        "Summary": {
            "rich_text": [
                {
                    "text": {
                        "content": (
                            f"Internal {omics_type.lower()} dataset ingested from {file_path}. "
                            f"Contains {raw_rows} raw entries, {entity_count} unique {entity_name}."
                        )
                    }
                }
            ]
        },
    }

    # Optional: Set Omics Type if property exists
    try:
        props["Omics Type"] = {"select": {"name": omics_type}}
    except Exception:
        # Property may not exist - that's okay
        pass

    # Create page
    payload = {
        "parent": {"database_id": cfg.exp_data_db_id},
        "properties": props,
    }

    try:
        resp = requests.post(
            f"{cfg.base_url}/pages",
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        page_id = resp.json()["id"]
        logger.info(
            "[INGEST][%s] Created new dataset page %s for file %s",
            omics_type.upper(),
            page_id,
            file_path,
        )
        return page_id
    except Exception as e:
        # Check if error is due to Omics Type property not existing
        if hasattr(e, "response") and e.response is not None:
            error_text = e.response.text
            if "Omics Type" in error_text and "not a property" in error_text:
                # Retry without Omics Type property
                logger.info(
                    "[INGEST][%s] Omics Type property not found, creating page without it",
                    omics_type.upper(),
                )
                if "Omics Type" in props:
                    del props["Omics Type"]
                payload_retry = {
                    "parent": {"database_id": cfg.exp_data_db_id},
                    "properties": props,
                }
                resp = requests.post(
                    f"{cfg.base_url}/pages",
                    headers=notion_headers(),
                    json=payload_retry,
                    timeout=30,
                )
                resp.raise_for_status()
                page_id = resp.json()["id"]
                logger.info(
                    "[INGEST][%s] Created new dataset page %s for file %s (without Omics Type)",
                    omics_type.upper(),
                    page_id,
                    file_path,
                )
                return page_id
            else:
                logger.error(
                    "[INGEST][%s] Error creating dataset page: %s - Response: %s",
                    omics_type.upper(),
                    str(e),
                    error_text,
                )
                raise
        else:
            logger.error(
                "[INGEST][%s] Error creating dataset page: %r",
                omics_type.upper(),
                e,
            )
            raise


def attach_file_to_page(
    page_id: str,
    file_path: str,
    omics_type: str,
) -> None:
    """
    Attach a file to a Notion page.

    Args:
        page_id: Notion page ID
        file_path: Path to file to attach
        omics_type: Type of omics (for logging prefix)
    """
    # Note: File attachment via Notion API requires uploading to external storage first
    # For now, we'll add a note in the Summary about the file location
    logger.info(
        "[INGEST][%s] File attachment for %s noted (file: %s)",
        omics_type.upper(),
        page_id,
        file_path,
    )


def link_to_programs_and_experiments(
    dataset_page_id: str,
    omics_type: str,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Link dataset to Programs and Experiments via Notion relations.

    Args:
        dataset_page_id: Dataset page ID
        omics_type: Type of omics (for logging prefix)
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    # Link to Programs
    if program_ids:
        for program_id in program_ids:
            try:
                logger.info(
                    "[INGEST][%s] Linking dataset %s to program %s",
                    omics_type.upper(),
                    dataset_page_id,
                    program_id,
                )
                # TODO: Implement actual relation updates when schema is confirmed
            except Exception as e:
                logger.warning(
                    "[INGEST][%s] Error linking to program %s: %r",
                    omics_type.upper(),
                    program_id,
                    e,
                )

    # Link to Experiments
    if experiment_ids:
        for experiment_id in experiment_ids:
            try:
                logger.info(
                    "[INGEST][%s] Linking dataset %s to experiment %s",
                    omics_type.upper(),
                    dataset_page_id,
                    experiment_id,
                )
                # TODO: Implement actual relation updates when schema is confirmed
            except Exception as e:
                logger.warning(
                    "[INGEST][%s] Error linking to experiment %s: %r",
                    omics_type.upper(),
                    experiment_id,
                    e,
                )

