"""
Batch signature scoring with caching support.

Provides efficient batch scoring of multiple datasets against signatures
using feature caching to minimize Notion API calls.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
from amprenta_rag.ingestion.multi_omics_scoring import (
    extract_dataset_features_by_type,
)
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset,
    SignatureMatchResult,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def score_datasets_against_signatures_batch(
    dataset_page_ids: List[str],
    signature_page_ids: Optional[List[str]] = None,
    overlap_threshold: float = 0.3,
    preload_cache: bool = True,
    use_cache: bool = True,
) -> Dict[str, List[SignatureMatchResult]]:
    """
    Score multiple datasets against signatures efficiently using caching.

    This function:
    1. Pre-loads all dataset features into cache (if preload_cache=True)
    2. Scores each dataset against signatures using cached features
    3. Returns results for all datasets

    Args:
        dataset_page_ids: List of dataset page IDs to score
        signature_page_ids: Optional list of signature page IDs (if None, uses all)
        overlap_threshold: Minimum overlap for matches (default: 0.3)
        preload_cache: Whether to preload dataset features into cache (default: True)
        use_cache: Whether to use cache during scoring (default: True)

    Returns:
        Dictionary mapping dataset_page_id → list of SignatureMatchResult objects
    """
    # Parameter retained for backward compatibility; filtering is not implemented in this module.
    del signature_page_ids
    results: Dict[str, List[SignatureMatchResult]] = {}

    logger.info(
        "[INGEST][BATCH-SCORE] Starting batch signature scoring for %d datasets",
        len(dataset_page_ids),
    )

    # Preload cache if requested
    if preload_cache and use_cache:
        logger.info(
            "[INGEST][BATCH-SCORE] Preloading features for %d datasets into cache",
            len(dataset_page_ids),
        )
        cache = get_feature_cache()
        preload_results = cache.preload_datasets(
            dataset_page_ids=dataset_page_ids,
            extract_fn=extract_dataset_features_by_type,
        )
        preload_success = sum(1 for v in preload_results.values() if v)
        logger.info(
            "[INGEST][BATCH-SCORE] Successfully preloaded %d/%d datasets",
            preload_success,
            len(dataset_page_ids),
        )

    # Score each dataset
    for idx, dataset_page_id in enumerate(dataset_page_ids, 1):
        try:
            logger.info(
                "[INGEST][BATCH-SCORE] Scoring dataset %d/%d: %s",
                idx,
                len(dataset_page_ids),
                dataset_page_id,
            )

            matches = find_matching_signatures_for_dataset(
                dataset_page_id=dataset_page_id,
                overlap_threshold=overlap_threshold,
            )

            results[dataset_page_id] = matches

            logger.info(
                "[INGEST][BATCH-SCORE] Dataset %s: Found %d matching signatures",
                dataset_page_id,
                len(matches),
            )

        except Exception as e:
            logger.warning(
                "[INGEST][BATCH-SCORE] Error scoring dataset %s: %r",
                dataset_page_id,
                e,
            )
            results[dataset_page_id] = []

    logger.info(
        "[INGEST][BATCH-SCORE] Batch scoring complete: %d datasets processed",
        len(dataset_page_ids),
    )

    return results


def clear_cache_for_datasets(dataset_page_ids: Optional[List[str]] = None) -> None:
    """
    Clear cache entries for specific datasets or all datasets.

    Args:
        dataset_page_ids: Optional list of dataset page IDs to clear.
                         If None, clears entire cache.
    """
    cache = get_feature_cache()

    if dataset_page_ids:
        for dataset_page_id in dataset_page_ids:
            cache.invalidate(dataset_page_id)
        logger.info(
            "[INGEST][BATCH-SCORE] Cleared cache for %d datasets",
            len(dataset_page_ids),
        )
    else:
        cache.clear()
        logger.info(
            "[INGEST][BATCH-SCORE] Cleared entire feature cache",
        )


def get_cache_stats() -> Dict[str, int]:
    """
    Get cache statistics.

    Returns:
        Dictionary with cache statistics
    """
    cache = get_feature_cache()
    return cache.get_stats()

