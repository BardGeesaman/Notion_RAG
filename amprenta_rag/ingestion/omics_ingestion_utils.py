"""
Shared utilities for omics dataset ingestion pipelines.

This module provides common functions used across lipidomics, metabolomics,
proteomics, and transcriptomics ingestion modules to reduce code duplication.

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

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
    Stub: Notion support removed. Returns empty string.
    
    Previously created a new Experimental Data Asset page in Notion.
    Datasets are now created directly in Postgres.
    """
    logger.debug(
        "[OMICS-INGESTION-UTILS] create_omics_dataset_page() is a no-op (Notion removed)"
    )
    return ""


def attach_file_to_page(
    page_id: str,
    file_path: str,
    omics_type: str,
) -> None:
    """
    Stub: Notion support removed. No-op.
    
    Previously attached a file to a Notion page.
    """
    logger.debug(
        "[OMICS-INGESTION-UTILS] attach_file_to_page() is a no-op (Notion removed)"
    )


def link_to_programs_and_experiments(
    dataset_page_id: str,
    omics_type: str,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Stub: Notion support removed. No-op.
    
    Previously linked dataset to Programs and Experiments via Notion relations.
    Linking is now done directly in Postgres.
    """
    logger.debug(
        "[OMICS-INGESTION-UTILS] link_to_programs_and_experiments() is a no-op (Notion removed)"
    )

