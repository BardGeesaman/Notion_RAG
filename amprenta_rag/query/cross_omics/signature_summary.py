"""
Cross-omics signature summary generation.

Generates multi-omics summaries for Signatures by aggregating evidence
from matched datasets and signature components.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Signature as SignatureModel
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.signature_summary_postgres import (
    cross_omics_signature_summary_postgres,
)

logger = get_logger(__name__)


def cross_omics_signature_summary(
    signature_page_id: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Generate a cross-omics summary for a Signature (Postgres-backed).
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for signature %s",
        signature_page_id,
    )

    # Try UUID first
    try:
        signature_uuid = UUID(signature_page_id)
        return cross_omics_signature_summary_postgres(
            signature_id=signature_uuid,
            top_k_datasets=top_k_datasets,
            top_k_chunks=top_k_chunks,
        )
    except ValueError:
        pass

    # Fallback: lookup by legacy notion_page_id
    db = next(get_db())
    try:
        signature: Optional[SignatureModel] = (
            db.query(SignatureModel)
            .filter(SignatureModel.notion_page_id == signature_page_id)
            .first()
        )
        if not signature:
            return f"Error: Signature with ID {signature_page_id} not found in Postgres."

        return cross_omics_signature_summary_postgres(
            signature_id=signature.id,
            top_k_datasets=top_k_datasets,
            top_k_chunks=top_k_chunks,
        )
    finally:
        db.close()

