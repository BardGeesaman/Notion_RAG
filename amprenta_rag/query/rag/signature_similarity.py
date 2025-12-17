"""
Signature similarity query helpers.

This module exists primarily for backward compatibility: older CLI/tools import
``signature_similarity_query`` from the RAG query engine, which historically
accepted a dataset "page id". In the Postgres-first architecture, datasets are
identified by UUIDs.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.ingestion.postgres_signature_matching import (
    find_matching_signatures_for_postgres_dataset,
)


def signature_similarity_query(
    *,
    dataset_page_id: str,
    top_k: int = 10,
    overlap_threshold: float = 0.3,
    omics_type: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Return top matching signatures for a dataset.

    Args:
        dataset_page_id: For compatibility, a string identifier. In Postgres mode
            this should be a UUID string of the dataset.
        top_k: Maximum results to return.
        overlap_threshold: Minimum overlap to include a signature match.
        omics_type: Optional modality hint.

    Returns:
        List of dicts compatible with the legacy CLI output.
    """
    dataset_id = UUID(dataset_page_id)
    matches = find_matching_signatures_for_postgres_dataset(
        dataset_id=dataset_id,
        overlap_threshold=overlap_threshold,
        omics_type=omics_type,
    )

    # Sort by score (desc) then overlap (desc)
    matches_sorted = sorted(
        matches,
        key=lambda m: (float(getattr(m, "score", 0.0)), float(getattr(m, "overlap_fraction", 0.0))),
        reverse=True,
    )[: max(0, int(top_k))]

    out: List[Dict[str, Any]] = []
    for m in matches_sorted:
        out.append(
            {
                "signature_page_id": getattr(m, "signature_page_id", ""),
                "signature_name": getattr(m, "signature_name", "Unknown"),
                "score": float(getattr(m, "score", 0.0)),
                "overlap_fraction": float(getattr(m, "overlap_fraction", 0.0)),
                "matched_components": sorted(list(getattr(m, "matched_components", []) or [])),
                "missing_components": sorted(list(getattr(m, "missing_components", []) or [])),
                "conflicting_components": sorted(list(getattr(m, "conflicting_components", []) or [])),
            }
        )
    return out


