"""
Pathway enrichment analysis utilities.

Provides statistical enrichment analysis functions for pathway analysis.
"""

from __future__ import annotations

from typing import List, Optional, Set

from amprenta_rag.analysis.pathway_analysis import (
    PathwayEnrichmentResult,
    perform_pathway_enrichment,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def enrich_dataset_pathways(
    dataset_page_id: str,
    p_value_threshold: float = 0.05,
    pathway_sources: Optional[List[str]] = None,
) -> List[PathwayEnrichmentResult]:
    """
    Perform pathway enrichment for a dataset.

    Args:
        dataset_page_id: Notion page ID of dataset
        p_value_threshold: P-value threshold for significance
        pathway_sources: List of sources to use (["KEGG", "Reactome"] or None for all)

    Returns:
        List of PathwayEnrichmentResult objects
    """
    logger.info(
        "[ANALYSIS][ENRICHMENT] Performing pathway enrichment for dataset %s",
        dataset_page_id,
    )

    # Extract features from dataset
    from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type

    features_by_type = extract_dataset_features_by_type(
        dataset_page_id,
        use_cache=True,
    )

    # Combine all features
    all_features: Set[str] = set()
    feature_types: Set[str] = set()

    for omics_type, features in features_by_type.items():
        all_features.update(features)
        if features:
            feature_types.add(omics_type)

    if not all_features:
        logger.warning(
            "[ANALYSIS][ENRICHMENT] No features found for dataset %s",
            dataset_page_id,
        )
        return []

    # Perform enrichment
    results = perform_pathway_enrichment(
        input_features=all_features,
        input_feature_types=feature_types,
        pathway_sources=pathway_sources,
        p_value_threshold=p_value_threshold,
    )

    logger.info(
        "[ANALYSIS][ENRICHMENT] Found %d enriched pathways for dataset %s",
        len(results),
        dataset_page_id,
    )

    return results


def enrich_signature_pathways(
    signature_page_id: str,
    p_value_threshold: float = 0.05,
    pathway_sources: Optional[List[str]] = None,
) -> List[PathwayEnrichmentResult]:
    """
    Perform pathway enrichment for a signature.

    Args:
        signature_page_id: Notion page ID of signature
        p_value_threshold: P-value threshold for significance
        pathway_sources: List of sources to use (["KEGG", "Reactome"] or None for all)

    Returns:
        List of PathwayEnrichmentResult objects
    """
    logger.info(
        "[ANALYSIS][ENRICHMENT] Performing pathway enrichment for signature %s",
        signature_page_id,
    )

    # Load signature
    from amprenta_rag.ingestion.signature_matching.signature_loader import (
        fetch_all_signatures_from_notion,
        load_signature_from_notion_page,
    )

    all_signature_pages = fetch_all_signatures_from_notion()
    signature_page = next(
        (p for p in all_signature_pages if p.get("id") == signature_page_id),
        None,
    )

    if not signature_page:
        logger.warning(
            "[ANALYSIS][ENRICHMENT] Signature %s not found",
            signature_page_id,
        )
        return []

    signature = load_signature_from_notion_page(signature_page)
    if not signature:
        logger.warning(
            "[ANALYSIS][ENRICHMENT] Could not load signature %s",
            signature_page_id,
        )
        return []

    # Extract features from signature
    all_features: Set[str] = set()
    feature_types: Set[str] = set()

    for comp in signature.components:
        all_features.add(comp.feature_name)
        if comp.feature_type:
            feature_types.add(comp.feature_type)

    if not all_features:
        logger.warning(
            "[ANALYSIS][ENRICHMENT] No features found in signature %s",
            signature_page_id,
        )
        return []

    # Perform enrichment
    results = perform_pathway_enrichment(
        input_features=all_features,
        input_feature_types=feature_types,
        pathway_sources=pathway_sources,
        p_value_threshold=p_value_threshold,
    )

    logger.info(
        "[ANALYSIS][ENRICHMENT] Found %d enriched pathways for signature %s",
        len(results),
        signature_page_id,
    )

    return results

