"""
Batch and async-optimized feature linking.

Provides high-performance feature linking using:
- Batch API calls for feature page lookups
- Batch relation updates
- Async/parallel processing
- Caching to avoid redundant API calls
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Global cache for feature page lookups (in-memory, per-session)
_feature_page_cache: Dict[Tuple[str, str], Optional[str]] = {}
_relation_property_cache: Dict[str, Optional[str]] = {}  # feature_page_id -> relation_property_name


def batch_find_or_create_feature_pages(
    features: List[Tuple[str, str]],
    max_workers: int = 10,
) -> Dict[Tuple[str, str], Optional[str]]:
    """
    Batch find or create feature pages.

    Args:
        features: List of (feature_type, feature_name) tuples
        max_workers: Maximum number of parallel workers

    Returns:
        Dictionary mapping (feature_type, feature_name) -> page_id (or None if failed)
    """
    logger.info(
        "[INGEST][FEATURE][BATCH] Batch processing %d features (max_workers=%d)",
        len(features),
        max_workers,
    )

    # Check cache first
    results: Dict[Tuple[str, str], Optional[str]] = {}
    uncached_features: List[Tuple[str, str]] = []

    for feature_tuple in features:
        if feature_tuple in _feature_page_cache:
            results[feature_tuple] = _feature_page_cache[feature_tuple]
        else:
            uncached_features.append(feature_tuple)

    if not uncached_features:
        logger.debug(
            "[INGEST][FEATURE][BATCH] All %d features found in cache",
            len(features),
        )
        return results

    logger.debug(
        "[INGEST][FEATURE][BATCH] %d features in cache, %d need processing",
        len(features) - len(uncached_features),
        len(uncached_features),
    )

    # Process uncached features in parallel
    from amprenta_rag.ingestion.features.general_linking import _find_or_create_feature_page

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_feature = {
            executor.submit(_find_or_create_feature_page, feat_type, feat_name): (feat_type, feat_name)
            for feat_type, feat_name in uncached_features
        }

        for future in as_completed(future_to_feature):
            feature_tuple = future_to_feature[future]
            try:
                page_id = future.result()
                results[feature_tuple] = page_id
                _feature_page_cache[feature_tuple] = page_id  # Cache result
            except Exception as e:
                logger.warning(
                    "[INGEST][FEATURE][BATCH] Error processing feature %s: %r",
                    feature_tuple,
                    e,
                )
                results[feature_tuple] = None

    logger.info(
        "[INGEST][FEATURE][BATCH] Batch processed %d features: %d succeeded, %d failed",
        len(uncached_features),
        sum(1 for v in results.values() if v),
        sum(1 for v in results.values() if not v),
    )

    return results


def batch_add_dataset_relations(
    feature_page_ids: List[Tuple[str, str, str]],
    dataset_page_id: str,
    max_workers: int = 10,
) -> None:
    """
    Batch add dataset relations to multiple feature pages.

    Args:
        feature_page_ids: List of (feature_type, feature_name, feature_page_id) tuples
        dataset_page_id: Dataset page ID to link
        max_workers: Maximum number of parallel workers
    """
    if not feature_page_ids:
        return

    logger.info(
        "[INGEST][FEATURE][BATCH] Batch adding dataset relations to %d features (max_workers=%d)",
        len(feature_page_ids),
        max_workers,
    )

    from amprenta_rag.ingestion.features.general_linking import _add_dataset_relation  # type: ignore[attr-defined]

    # Process in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for feature_type, feature_name, page_id in feature_page_ids:
            if page_id:  # Only process if page was found/created
                future = executor.submit(
                    _add_dataset_relation,
                    page_id,
                    dataset_page_id,
                    feature_type,
                )
                futures.append(future)

        # Wait for all to complete (errors are logged but don't stop processing)
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.debug(
                    "[INGEST][FEATURE][BATCH] Error in batch relation update: %r",
                    e,
                )
                # Continue - individual errors are non-blocking

    logger.info(
        "[INGEST][FEATURE][BATCH] Completed batch relation updates for %d features",
        len(feature_page_ids),
    )


def batch_link_features(
    features: List[Tuple[str, str]],
    dataset_page_id: str,
    max_workers: int = 10,
    enable_linking: bool = True,
) -> Dict[Tuple[str, str], Optional[str]]:
    """
    Batch link multiple features to a dataset.

    This is the high-performance version of link_feature() that processes
    multiple features in parallel with batching.

    Args:
        features: List of (feature_type, feature_name) tuples
        dataset_page_id: Dataset page ID to link features to
        max_workers: Maximum number of parallel workers
        enable_linking: If False, skip linking (only find/create pages)

    Returns:
        Dictionary mapping (feature_type, feature_name) -> page_id
    """
    logger.info(
        "[INGEST][FEATURE][BATCH] Batch linking %d features to dataset %s",
        len(features),
        dataset_page_id,
    )

    # Step 1: Batch find/create feature pages
    feature_pages = batch_find_or_create_feature_pages(features, max_workers=max_workers)

    # Step 2: Batch add dataset relations
    if enable_linking:
        feature_page_ids = [
            (feat_type, feat_name, page_id) for (feat_type, feat_name), page_id in feature_pages.items() if page_id
        ]
        batch_add_dataset_relations(feature_page_ids, dataset_page_id, max_workers=max_workers)

    return feature_pages


def clear_feature_cache() -> None:
    """Clear the in-memory feature page cache."""
    global _feature_page_cache, _relation_property_cache
    _feature_page_cache.clear()
    _relation_property_cache.clear()
    logger.debug("[INGEST][FEATURE][BATCH] Cleared feature cache")
