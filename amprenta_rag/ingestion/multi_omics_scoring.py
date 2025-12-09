"""
Multi-omics signature scoring engine.

Extends signature scoring to work across all omics types:
- Gene signatures → Transcriptomics datasets
- Protein signatures → Proteomics datasets  
- Metabolite signatures → Metabolomics datasets
- Lipid signatures → Lipidomics datasets
- Mixed signatures → Multiple dataset types
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset, Feature, dataset_feature_assoc
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature
from amprenta_rag.signatures.signature_scoring import SignatureScoreResult

logger = get_logger(__name__)

def extract_dataset_features_by_type(
    dataset_page_id: str,
    omics_type: Optional[str] = None,
    use_cache: bool = True,
    force_refresh: bool = False,
) -> Dict[str, Set[str]]:
    """
    Extract features from Postgres grouped by feature type.

    Args:
        dataset_page_id: Dataset identifier (UUID or legacy Notion page ID)
        omics_type: Optional omics type hint (unused in Postgres path)
        use_cache: Whether to use feature cache (default: True)
        force_refresh: Force refresh even if cached (default: False)

    Returns:
        Dictionary mapping feature_type → set of feature names.
    """
    if use_cache and not force_refresh:
        from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache

        cache = get_feature_cache()
        cached_features = cache.get_features(dataset_page_id)
        if cached_features:
            logger.debug(
                "[INGEST][MULTI-OMICS-SCORE] Using cached features for dataset %s",
                dataset_page_id,
            )
            return cached_features

    features_by_type: Dict[str, Set[str]] = {
        "gene": set(),
        "protein": set(),
        "metabolite": set(),
        "lipid": set(),
    }

    db_gen = get_db()
    db: Session = next(db_gen)
    try:
        dataset: Optional[Dataset] = None

        # Try UUID lookup first
        try:
            dataset_uuid = UUID(dataset_page_id)
            dataset = db.query(Dataset).filter(Dataset.id == dataset_uuid).first()
        except ValueError:
            dataset = None

        if dataset is None:
            dataset = (
                db.query(Dataset)
                .filter(Dataset.notion_page_id == dataset_page_id)
                .first()
            )

        if dataset is None:
            logger.warning(
                "[INGEST][MULTI-OMICS-SCORE] Dataset not found for id=%s",
                dataset_page_id,
            )
            return features_by_type

        rows = (
            db.query(Feature.feature_type, Feature.name)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id == dataset.id)
            .all()
        )

        for feature_type, name in rows:
            ft = (feature_type or "").lower()
            if ft in features_by_type:
                features_by_type[ft].add(name)

        total_features = sum(len(v) for v in features_by_type.values())
        logger.info(
            "[INGEST][MULTI-OMICS-SCORE] Extracted %d total features for dataset %s",
            total_features,
            dataset_page_id,
        )
    finally:
        db.close()

    if use_cache:
        from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache

        cache = get_feature_cache()
        cache.set_features(
            dataset_id=dataset_page_id,
            features=features_by_type,
        )

    return features_by_type


def score_multi_omics_signature_against_dataset(
    signature: Signature,
    dataset_features_by_type: Dict[str, Set[str]],
    dataset_directions: Optional[Dict[str, Dict[str, str]]] = None,
) -> SignatureScoreResult:
    """
    Score a multi-omics signature against dataset features.

    Matches signature components to dataset features based on feature_type.

    Args:
        signature: Multi-omics signature definition
        dataset_features_by_type: Dict mapping feature_type → set of feature names
        dataset_directions: Optional dict mapping feature_type → {feature_name → direction}

    Returns:
        SignatureScoreResult with scoring details
    """
    from amprenta_rag.signatures.signature_scoring import ComponentMatch

    if dataset_directions is None:
        dataset_directions = {}

    component_matches: List[ComponentMatch] = []
    total_weighted_score = 0.0
    total_weight = 0.0

    # Import normalization functions
    from amprenta_rag.ingestion.transcriptomics_ingestion import normalize_gene_identifier
    from amprenta_rag.ingestion.proteomics_ingestion import normalize_protein_identifier
    from amprenta_rag.ingestion.metabolomics_ingestion import normalize_metabolite_name
    from amprenta_rag.ingestion.lipidomics_ingestion import normalize_lipid_species
    from amprenta_rag.signatures.species_matching import normalize_species_name

    for comp in signature.components:
        feature_type = getattr(comp, "feature_type", "lipid")  # Default to lipid
        feature_name = getattr(comp, "feature_name", comp.species)

        # Get dataset features of matching type
        dataset_features = dataset_features_by_type.get(feature_type, set())

        # Normalize feature name based on type
        normalized_sig_name = feature_name
        if feature_type == "gene":
            normalized_sig_name = normalize_gene_identifier(feature_name)
        elif feature_type == "protein":
            normalized_sig_name = normalize_protein_identifier(feature_name)
        elif feature_type == "metabolite":
            normalized_sig_name = normalize_metabolite_name(feature_name)
        elif feature_type == "lipid":
            normalized = normalize_lipid_species(feature_name)
            if normalized:
                normalized_sig_name = normalized
            else:
                normalized_sig_name = normalize_species_name(feature_name)

        # Normalize all dataset features for comparison
        normalized_dataset_features = {}
        for dataset_feature in dataset_features:
            if feature_type == "gene":
                normalized = normalize_gene_identifier(dataset_feature)
            elif feature_type == "protein":
                normalized = normalize_protein_identifier(dataset_feature)
            elif feature_type == "metabolite":
                normalized = normalize_metabolite_name(dataset_feature)
            elif feature_type == "lipid":
                normalized_val = normalize_lipid_species(dataset_feature)
                if normalized_val:
                    normalized = normalized_val
                else:
                    normalized = normalize_species_name(dataset_feature)
            else:
                normalized = dataset_feature

            normalized_dataset_features[normalized] = dataset_feature

        # Try to find match
        matched_dataset_feature = None
        match_type = "none"

        # Exact match
        if normalized_sig_name in normalized_dataset_features:
            matched_dataset_feature = normalized_dataset_features[normalized_sig_name]
            match_type = "exact"
        else:
            # Try case-insensitive match
            normalized_sig_lower = normalized_sig_name.lower()
            for norm_name, orig_name in normalized_dataset_features.items():
                if norm_name.lower() == normalized_sig_lower:
                    matched_dataset_feature = orig_name
                    match_type = "exact"
                    break

        # Determine direction match
        direction_match = "unknown"
        match_value = 0.0

        if matched_dataset_feature:
            # Get dataset direction (if available)
            type_directions = dataset_directions.get(feature_type, {})
            dataset_dir = type_directions.get(matched_dataset_feature)
            sig_dir = comp.direction

            if sig_dir and dataset_dir:
                # Both have direction info
                if sig_dir == dataset_dir:
                    direction_match = "match"
                    match_value = 1.0
                elif (sig_dir == "↑" and dataset_dir == "↓") or (
                    sig_dir == "↓" and dataset_dir == "↑"
                ):
                    direction_match = "conflict"
                    match_value = -1.0
                else:
                    direction_match = "neutral"
                    match_value = 0.3
            elif sig_dir:
                # Only signature has direction (dataset present but no direction)
                direction_match = "neutral"
                match_value = 0.3
            else:
                # Species present but no direction info
                direction_match = "unknown"
                match_value = 0.3
        else:
            # Feature not found in dataset
            direction_match = "missing"
            match_value = 0.0

        # Weight the match
        weight = comp.weight or 1.0
        weighted_score = weight * match_value
        total_weighted_score += weighted_score
        total_weight += weight

        component_matches.append(
            ComponentMatch(
                signature_species=feature_name,
                matched_dataset_species=matched_dataset_feature,
                match_type=match_type,
                direction_match=direction_match,
                weight=weight,
            )
        )

        logger.debug(
            "[INGEST][MULTI-OMICS-SCORE] Component %s (%s): %s (match: %s, score: %.2f)",
            feature_name,
            feature_type,
            direction_match,
            "yes" if matched_dataset_feature else "no",
            match_value,
        )

    # Calculate final score
    if total_weight > 0:
        total_score = total_weighted_score / total_weight
        # Normalize to 0-1 range (since match_value can be -1 to +1)
        total_score = (total_score + 1.0) / 2.0
    else:
        total_score = 0.0

    # Categorize results
    missing_species = [
        comp.signature_species
        for comp in component_matches
        if comp.matched_dataset_species is None
    ]

    conflicting_species = [
        comp.signature_species
        for comp in component_matches
        if comp.direction_match == "conflict"
    ]

    matched_species = [
        comp.matched_dataset_species
        for comp in component_matches
        if comp.matched_dataset_species is not None
    ]

    logger.info(
        "[INGEST][MULTI-OMICS-SCORE] Scored signature '%s': score=%.3f, "
        "matched=%d/%d, missing=%d, conflicts=%d",
        signature.name,
        total_score,
        len(matched_species),
        len(signature.components),
        len(missing_species),
        len(conflicting_species),
    )

    return SignatureScoreResult(
        total_score=total_score,
        component_matches=component_matches,
        missing_species=missing_species,
        conflicting_species=conflicting_species,
        matched_species=matched_species,
    )

