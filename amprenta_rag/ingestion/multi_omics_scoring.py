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

from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature, SignatureComponent
from amprenta_rag.signatures.signature_scoring import SignatureScoreResult

logger = get_logger(__name__)


def extract_dataset_features_by_type(
    dataset_page_id: str,
    omics_type: Optional[str] = None,
    use_cache: bool = True,
    force_refresh: bool = False,
) -> Dict[str, Set[str]]:
    """
    Extract features from a dataset grouped by feature type.

    Queries feature databases for features that are linked to this dataset.
    Uses caching to avoid repeated Notion API calls.

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        omics_type: Optional omics type hint (Lipidomics, Metabolomics, Proteomics, Transcriptomics)
        use_cache: Whether to use feature cache (default: True)
        force_refresh: Force refresh even if cached (default: False)

    Returns:
        Dictionary mapping feature_type → set of feature names:
        {
            "gene": {"TP53", "TNF", ...},
            "protein": {"P04637", ...},
            "metabolite": {"Glutamate", ...},
            "lipid": {"Cer(d18:1/16:0)", ...}
        }
    """
    from amprenta_rag.clients.notion_client import notion_headers
    from amprenta_rag.config import get_config

    # Check cache first
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

    cfg = get_config()
    features_by_type: Dict[str, Set[str]] = {
        "gene": set(),
        "protein": set(),
        "metabolite": set(),
        "lipid": set(),
    }

    # Map feature types to database IDs and relation property names
    feature_db_map = {
        "gene": (
            cfg.notion.gene_features_db_id if hasattr(cfg.notion, "gene_features_db_id") else None,
            ["Transcriptomics Datasets", "Datasets", "Related Datasets"],
        ),
        "protein": (
            cfg.notion.protein_features_db_id if hasattr(cfg.notion, "protein_features_db_id") else None,
            ["Proteomics Datasets", "Datasets", "Related Datasets"],
        ),
        "metabolite": (
            cfg.notion.metabolite_features_db_id if hasattr(cfg.notion, "metabolite_features_db_id") else None,
            ["Metabolomics Datasets", "Datasets", "Related Datasets"],
        ),
        "lipid": (
            cfg.notion.lipid_species_db_id if hasattr(cfg.notion, "lipid_species_db_id") else None,
            ["Experimental Data Assets", "Datasets", "Related Datasets"],
        ),
    }

    # Query each feature database for features linked to this dataset
    for feature_type, (db_id, relation_property_candidates) in feature_db_map.items():
        if not db_id:
            logger.debug(
                "[INGEST][MULTI-OMICS-SCORE] %s database ID not configured, skipping",
                feature_type,
            )
            continue

        try:
            # Try each candidate relation property name
            found_features = False
            for relation_prop in relation_property_candidates:
                try:
                    url = f"{cfg.notion.base_url}/databases/{db_id}/query"
                    all_features = []
                    has_more = True
                    start_cursor = None

                    while has_more:
                        payload = {
                            "filter": {
                                "property": relation_prop,
                                "relation": {"contains": dataset_page_id},
                            },
                            "page_size": 100,
                        }
                        if start_cursor:
                            payload["start_cursor"] = start_cursor

                        resp = requests.post(
                            url,
                            headers=notion_headers(),
                            json=payload,
                            timeout=30,
                        )
                        resp.raise_for_status()

                        data = resp.json()
                        all_features.extend(data.get("results", []))
                        has_more = data.get("has_more", False)
                        start_cursor = data.get("next_cursor")

                    # Extract feature names from pages
                    if all_features:
                        for feature_page in all_features:
                            props = feature_page.get("properties", {}) or {}
                            
                            # Get feature name (usually "Name" title property)
                            name_prop = props.get("Name", {}).get("title", []) or []
                            if name_prop:
                                feature_name = name_prop[0].get("plain_text", "").strip()
                                if feature_name:
                                    features_by_type[feature_type].add(feature_name)
                                    found_features = True

                        if found_features:
                            logger.debug(
                                "[INGEST][MULTI-OMICS-SCORE] Found %d %s features for dataset %s (via property '%s')",
                                len(features_by_type[feature_type]),
                                feature_type,
                                dataset_page_id,
                                relation_prop,
                            )
                            break  # Found features with this property, no need to try others

                except Exception as e:
                    # Property doesn't exist or query failed, try next candidate
                    logger.debug(
                        "[INGEST][MULTI-OMICS-SCORE] Could not query %s via property '%s': %r",
                        feature_type,
                        relation_prop,
                        e,
                    )
                    continue

            # If no features found with standard properties, try dynamic search
            if not found_features:
                try:
                    url = f"{cfg.notion.base_url}/databases/{db_id}/query"
                    # Query all pages and check properties dynamically
                    all_pages = []
                    has_more = True
                    start_cursor = None

                    while has_more:
                        payload = {"page_size": 100}
                        if start_cursor:
                            payload["start_cursor"] = start_cursor

                        resp = requests.post(
                            url,
                            headers=notion_headers(),
                            json=payload,
                            timeout=30,
                        )
                        resp.raise_for_status()

                        data = resp.json()
                        all_pages.extend(data.get("results", []))
                        has_more = data.get("has_more", False)
                        start_cursor = data.get("next_cursor")

                    # Check each page for dataset relation
                    for page in all_pages:
                        props = page.get("properties", {}) or {}
                        for prop_name, prop_data in props.items():
                            if prop_data.get("type") == "relation":
                                relations = prop_data.get("relation", []) or []
                                if any(r.get("id") == dataset_page_id for r in relations):
                                    # This feature is linked to the dataset
                                    name_prop = props.get("Name", {}).get("title", []) or []
                                    if name_prop:
                                        feature_name = name_prop[0].get("plain_text", "").strip()
                                        if feature_name:
                                            features_by_type[feature_type].add(feature_name)
                                            found_features = True
                                            break

                    if found_features:
                        logger.debug(
                            "[INGEST][MULTI-OMICS-SCORE] Found %d %s features for dataset %s (via dynamic search)",
                            len(features_by_type[feature_type]),
                            feature_type,
                            dataset_page_id,
                        )

                except Exception as e:
                    logger.debug(
                        "[INGEST][MULTI-OMICS-SCORE] Error in dynamic search for %s features: %r",
                        feature_type,
                        e,
                    )

        except Exception as e:
            logger.warning(
                "[INGEST][MULTI-OMICS-SCORE] Error extracting %s features for dataset %s: %r",
                feature_type,
                dataset_page_id,
                e,
            )
            continue

    total_features = sum(len(features) for features in features_by_type.values())
    logger.info(
        "[INGEST][MULTI-OMICS-SCORE] Extracted %d total features for dataset %s: "
        "genes=%d, proteins=%d, metabolites=%d, lipids=%d",
        total_features,
        dataset_page_id,
        len(features_by_type["gene"]),
        len(features_by_type["protein"]),
        len(features_by_type["metabolite"]),
        len(features_by_type["lipid"]),
    )

    # Cache the extracted features
    if use_cache:
        from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache

        cache = get_feature_cache()
        cache.set_features(
            dataset_page_id=dataset_page_id,
            features_by_type=features_by_type,
            omics_type=omics_type,
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

