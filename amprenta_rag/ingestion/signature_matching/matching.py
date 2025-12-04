"""
Signature matching and scoring logic.

Handles finding matching signatures for datasets and scoring them.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

from amprenta_rag.ingestion.multi_omics_scoring import (
    extract_dataset_features_by_type,
    score_multi_omics_signature_against_dataset,
)
from amprenta_rag.ingestion.signature_matching.models import SignatureMatchResult
from amprenta_rag.ingestion.signature_matching.signature_loader import (
    fetch_all_signatures_from_notion,
    load_signature_from_notion_page,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature
from amprenta_rag.signatures.signature_scoring import (
    SignatureScoreResult,
    score_signature,
)

logger = get_logger(__name__)


def score_signature_against_dataset(
    signature: Signature,
    dataset_species: Set[str],
    dataset_directions: Optional[Dict[str, str]] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> SignatureScoreResult:
    """
    Score a signature against a dataset.

    Computes overlap/similarity score between a lipid signature and a dataset.

    Args:
        signature: Signature definition
        dataset_species: Set of species names found in the dataset
        dataset_directions: Optional dict mapping species → direction (↑/↓)
        metadata: Optional dataset metadata (for future use)

    Returns:
        SignatureScoreResult with detailed scoring information
    """
    return score_signature(
        signature=signature,
        dataset_species=dataset_species,
        dataset_directions=dataset_directions,
    )


def find_matching_signatures_for_dataset(
    dataset_species: Optional[Set[str]] = None,
    dataset_directions: Optional[Dict[str, str]] = None,
    overlap_threshold: float = 0.3,
    dataset_page_id: Optional[str] = None,
    omics_type: Optional[str] = None,
) -> List[SignatureMatchResult]:
    """
    Find all signatures that match a dataset above the overlap threshold.

    Supports both legacy (lipid-only) and multi-omics signatures.

    Args:
        dataset_species: Set of species names in the dataset (legacy, for backward compat)
        dataset_directions: Optional dict mapping species → direction
        overlap_threshold: Minimum overlap fraction to consider a match (default: 0.3)
        dataset_page_id: Optional Notion page ID of dataset (for multi-omics extraction)
        omics_type: Optional omics type hint (Lipidomics, Metabolomics, etc.)

    Returns:
        List of SignatureMatchResult objects for matching signatures
    """
    matches: List[SignatureMatchResult] = []

    # Extract dataset features if dataset_page_id is provided (multi-omics mode)
    dataset_features_by_type: Optional[Dict[str, Set[str]]] = None
    if dataset_page_id:
        try:
            # Use cache by default for performance
            dataset_features_by_type = extract_dataset_features_by_type(
                dataset_page_id=dataset_page_id,
                omics_type=omics_type,
                use_cache=True,
            )
            logger.info(
                "[INGEST][SIGNATURE-MATCH] Extracted features from dataset %s for multi-omics scoring",
                dataset_page_id,
            )
        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Error extracting features from dataset %s: %r. "
                "Falling back to legacy scoring.",
                dataset_page_id,
                e,
            )
            dataset_features_by_type = None

    # Fetch all signatures from Notion
    signature_pages = fetch_all_signatures_from_notion()

    if not signature_pages:
        logger.debug(
            "[INGEST][SIGNATURE-MATCH] No signatures found in Notion",
        )
        return matches

    for sig_page in signature_pages:
        try:
            # Load signature from Notion page
            signature = load_signature_from_notion_page(sig_page)
            if not signature:
                continue

            # Determine if this is a multi-omics signature
            is_multi_omics = (
                signature.modalities
                and len(signature.modalities) > 0
                and (
                    len(signature.modalities) > 1
                    or signature.modalities[0].lower() != "lipid"
                )
            ) or any(
                getattr(comp, "feature_type", "lipid") != "lipid"
                for comp in signature.components
            )

            # Use multi-omics scoring if we have feature extraction and it's a multi-omics signature
            if dataset_features_by_type and is_multi_omics:
                # Multi-omics scoring
                score_result = score_multi_omics_signature_against_dataset(
                    signature=signature,
                    dataset_features_by_type=dataset_features_by_type,
                    dataset_directions=None,  # TODO: Support directions by feature type
                )
            elif dataset_species:
                # Legacy scoring (lipid-only)
                score_result = score_signature_against_dataset(
                    signature=signature,
                    dataset_species=dataset_species,
                    dataset_directions=dataset_directions,
                )
            else:
                # Can't score without dataset features
                logger.debug(
                    "[INGEST][SIGNATURE-MATCH] Skipping signature %s: no dataset features available",
                    signature.name,
                )
                continue

            # Calculate overlap fraction
            total_components = len(signature.components)
            matched_count = len(score_result.matched_species)
            overlap_fraction = (
                matched_count / total_components if total_components > 0 else 0.0
            )

            # Check if above threshold
            if overlap_fraction >= overlap_threshold:
                props = sig_page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                signature_name = (
                    name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                )

                match_result = SignatureMatchResult(
                    signature_page_id=sig_page.get("id", ""),
                    signature_name=signature_name,
                    score=score_result.total_score,
                    overlap_fraction=overlap_fraction,
                    matched_components=score_result.matched_species,
                    missing_components=score_result.missing_species,
                    conflicting_components=score_result.conflicting_species,
                    score_result=score_result,
                )
                matches.append(match_result)

                logger.info(
                    "[INGEST][SIGNATURE-MATCH] Found match: %s (overlap: %.2f, score: %.3f)",
                    signature_name,
                    overlap_fraction,
                    score_result.total_score,
                )

        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Error processing signature %s: %r",
                sig_page.get("id", ""),
                e,
            )
            continue

    return matches

