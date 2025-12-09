"""
Cross-omics dataset summary generation.

Generates multi-omics summaries for Datasets by aggregating evidence
from linked signatures, experiments, and programs.
"""

from __future__ import annotations

from typing import Optional
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset as DatasetModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.dataset_summary_postgres import (
    cross_omics_dataset_summary_postgres,
)

logger = get_logger(__name__)


def cross_omics_dataset_summary(
    dataset_page_id: str,
    top_k_chunks: int = 100,
) -> str:
    """
    Summarize cross-omics context for a single dataset.

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for dataset %s",
        dataset_page_id,
    )

    # Try UUID first
    try:
        dataset_uuid = UUID(dataset_page_id)
        return cross_omics_dataset_summary_postgres(dataset_uuid, top_k_chunks=top_k_chunks)
    except ValueError:
        pass

    # Fallback: lookup by legacy notion_page_id
    db = next(get_db())
    try:
        dataset: Optional[DatasetModel] = (
            db.query(DatasetModel)
            .filter(DatasetModel.notion_page_id == dataset_page_id)
            .first()
        )
        if not dataset:
            return f"Error: Dataset with ID {dataset_page_id} not found in Postgres."

        return cross_omics_dataset_summary_postgres(dataset.id, top_k_chunks=top_k_chunks)
    finally:
        db.close()

