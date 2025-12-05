"""
Postgres-based signature matching.

Finds matching signatures for Postgres datasets without requiring Notion page IDs.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set
from uuid import UUID

from amprenta_rag.ingestion.multi_omics_scoring import (
    score_multi_omics_signature_against_dataset,
)
from amprenta_rag.ingestion.postgres_feature_extraction import (
    extract_dataset_features_by_type_from_postgres,
)
from amprenta_rag.ingestion.postgres_signature_loader import (
    fetch_all_signatures_from_postgres,
    load_signature_from_postgres,
)
from amprenta_rag.ingestion.signature_matching.models import SignatureMatchResult
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature
from amprenta_rag.signatures.signature_scoring import (
    SignatureScoreResult,
    score_signature,
)

logger = get_logger(__name__)


def find_matching_signatures_for_postgres_dataset(
    dataset_id: UUID,
    dataset_species: Optional[Set[str]] = None,
    dataset_directions: Optional[Dict[str, str]] = None,
    overlap_threshold: float = 0.3,
    omics_type: Optional[str] = None,
) -> List[SignatureMatchResult]:
    """
    Find all signatures that match a Postgres dataset above the overlap threshold.
    
    Uses Postgres dataset_id to extract features and match against signatures.
    Supports both legacy (lipid-only) and multi-omics signatures.
    
    Args:
        dataset_id: Postgres UUID of the dataset
        dataset_species: Set of species names (legacy, for backward compat)
        dataset_directions: Optional dict mapping species → direction
        overlap_threshold: Minimum overlap fraction to consider a match (default: 0.3)
        omics_type: Optional omics type hint (Lipidomics, Metabolomics, etc.)
        
    Returns:
        List of SignatureMatchResult objects for matching signatures
    """
    matches: List[SignatureMatchResult] = []
    
    # Extract dataset features from Postgres
    dataset_features_by_type: Optional[Dict[str, Set[str]]] = None
    try:
        dataset_features_by_type = extract_dataset_features_by_type_from_postgres(
            dataset_id=dataset_id,
            omics_type=omics_type,
        )
        
        if dataset_features_by_type:
            total_features = sum(len(s) for s in dataset_features_by_type.values())
            logger.info(
                "[POSTGRES-SIGNATURE-MATCH] Extracted %d features from dataset %s for multi-omics scoring",
                total_features,
                dataset_id,
            )
        else:
            logger.warning(
                "[POSTGRES-SIGNATURE-MATCH] No features extracted from dataset %s",
                dataset_id,
            )
    except Exception as e:
        logger.warning(
            "[POSTGRES-SIGNATURE-MATCH] Error extracting features from dataset %s: %r. "
            "Falling back to legacy scoring if dataset_species provided.",
            dataset_id,
            e,
        )
        dataset_features_by_type = None
    
    # Fetch all signatures from Postgres
    signature_models = fetch_all_signatures_from_postgres()
    
    if not signature_models:
        logger.debug(
            "[POSTGRES-SIGNATURE-MATCH] No signatures found in Postgres",
        )
        return matches
    
    logger.info(
        "[POSTGRES-SIGNATURE-MATCH] Matching dataset %s against %d signature(s)",
        dataset_id,
        len(signature_models),
    )
    
    for sig_model in signature_models:
        try:
            # Load signature from Postgres model
            signature = load_signature_from_postgres(sig_model)
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
                    "[POSTGRES-SIGNATURE-MATCH] Skipping signature %s: no dataset features available",
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
                # Use signature UUID as identifier (store as string for compatibility)
                signature_id_str = str(sig_model.id)
                
                match_result = SignatureMatchResult(
                    signature_page_id=signature_id_str,  # UUID string for compatibility
                    signature_name=sig_model.name,
                    score=score_result.total_score,
                    overlap_fraction=overlap_fraction,
                    matched_components=score_result.matched_species,
                    missing_components=score_result.missing_species,
                    conflicting_components=score_result.conflicting_species,
                    score_result=score_result,
                )
                matches.append(match_result)
                
                logger.info(
                    "[POSTGRES-SIGNATURE-MATCH] Found match: %s (overlap: %.2f, score: %.3f)",
                    sig_model.name,
                    overlap_fraction,
                    score_result.total_score,
                )
        
        except Exception as e:
            logger.warning(
                "[POSTGRES-SIGNATURE-MATCH] Error processing signature %s: %r",
                sig_model.id,
                e,
            )
            continue
    
    logger.info(
        "[POSTGRES-SIGNATURE-MATCH] Found %d matching signature(s) for dataset %s",
        len(matches),
        dataset_id,
    )
    
    return matches

