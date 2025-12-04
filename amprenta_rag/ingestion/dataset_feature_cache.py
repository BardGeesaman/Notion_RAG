"""
Dataset feature cache for performance optimization.

This module provides a TTL-based caching system for dataset features extracted from Notion.
The cache significantly improves performance for batch signature scoring operations by
reducing redundant API calls.

Key Classes:
    - DatasetFeatureCache: Main cache class with TTL-based expiration

Key Functions:
    - get_feature_cache: Returns singleton cache instance
    - clear_feature_cache: Clears cache entries

Example:
    >>> from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
    >>> cache = get_feature_cache(ttl_seconds=3600)
    >>> features = cache.get_features("dataset-page-id-123")
    >>> if features is None:
    ...     features = extract_dataset_features_by_type("dataset-page-id-123")
    ...     cache.set_features("dataset-page-id-123", features)
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any, Callable, Dict, List, Optional, Set

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class DatasetFeatureCache:
    """
    In-memory cache for dataset feature sets to avoid repeated Notion API calls.

    Cache structure:
    {
        dataset_page_id: {
            "features_by_type": {
                "gene": set(...),
                "protein": set(...),
                "metabolite": set(...),
                "lipid": set(...)
            },
            "directions": {
                "gene": {"TP53": "↑", ...},
                ...
            },
            "timestamp": datetime,
            "omics_type": str
        }
    }
    """

    def __init__(self, ttl_seconds: int = 3600):
        """
        Initialize the cache.

        Args:
            ttl_seconds: Time-to-live for cache entries in seconds (default: 1 hour)
        """
        self.cache: Dict[str, Dict[str, Any]] = {}
        self.ttl = timedelta(seconds=ttl_seconds)
        logger.debug(
            "[CACHE] Initialized DatasetFeatureCache with TTL=%d seconds",
            ttl_seconds,
        )

    def get_features(
        self,
        dataset_page_id: str,
        force_refresh: bool = False,
    ) -> Optional[Dict[str, Set[str]]]:
        """
        Get cached features for a dataset, or None if not cached/stale.

        Args:
            dataset_page_id: Notion page ID of dataset (with dashes)
            force_refresh: If True, ignore cache and return None

        Returns:
            Dictionary mapping feature_type → set of feature names, or None if not cached/stale
        """
        if force_refresh:
            logger.debug(
                "[CACHE] Force refresh requested for dataset %s",
                dataset_page_id,
            )
            return None

        if dataset_page_id not in self.cache:
            logger.debug(
                "[CACHE] No cache entry for dataset %s",
                dataset_page_id,
            )
            return None

        entry = self.cache[dataset_page_id]
        age = datetime.now() - entry["timestamp"]

        if age > self.ttl:
            logger.debug(
                "[CACHE] Cache entry for dataset %s is stale (age: %s), will refresh",
                dataset_page_id,
                age,
            )
            return None

        logger.debug(
            "[CACHE] Cache hit for dataset %s (age: %s)",
            dataset_page_id,
            age,
        )
        return entry["features_by_type"]

    def get_directions(
        self,
        dataset_page_id: str,
    ) -> Optional[Dict[str, Dict[str, str]]]:
        """
        Get cached directions for a dataset.

        Args:
            dataset_page_id: Notion page ID of dataset (with dashes)

        Returns:
            Dictionary mapping feature_type → {feature_name: direction}, or None
        """
        if dataset_page_id not in self.cache:
            return None

        entry = self.cache[dataset_page_id]
        age = datetime.now() - entry["timestamp"]

        if age > self.ttl:
            return None

        return entry.get("directions")

    def set_features(
        self,
        dataset_page_id: str,
        features_by_type: Dict[str, Set[str]],
        directions: Optional[Dict[str, Dict[str, str]]] = None,
        omics_type: Optional[str] = None,
    ) -> None:
        """
        Store features in cache.

        Args:
            dataset_page_id: Notion page ID of dataset (with dashes)
            features_by_type: Dictionary mapping feature_type → set of feature names
            directions: Optional dictionary mapping feature_type → {feature_name: direction}
            omics_type: Optional omics type hint
        """
        feature_type_count = len([k for k, v in features_by_type.items() if v])
        total_features = sum(len(v) for v in features_by_type.values())

        self.cache[dataset_page_id] = {
            "features_by_type": features_by_type,
            "directions": directions or {},
            "timestamp": datetime.now(),
            "omics_type": omics_type,
        }

        logger.debug(
            "[CACHE] Cached features for dataset %s (%d feature types, %d total features)",
            dataset_page_id,
            feature_type_count,
            total_features,
        )

    def clear(self, dataset_page_id: Optional[str] = None) -> None:
        """
        Clear cache entry(ies).

        Args:
            dataset_page_id: Optional specific dataset to clear. If None, clears entire cache.
        """
        if dataset_page_id:
            if dataset_page_id in self.cache:
                del self.cache[dataset_page_id]
                logger.debug(
                    "[CACHE] Cleared cache for dataset %s",
                    dataset_page_id,
                )
            else:
                logger.debug(
                    "[CACHE] Dataset %s not in cache, nothing to clear",
                    dataset_page_id,
                )
        else:
            count = len(self.cache)
            self.cache.clear()
            logger.debug(
                "[CACHE] Cleared entire cache (%d entries)",
                count,
            )

    def get_stats(self) -> Dict[str, int]:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        total_features = sum(
            len(features)
            for entry in self.cache.values()
            for features in entry["features_by_type"].values()
        )

        return {
            "cached_datasets": len(self.cache),
            "total_features": total_features,
        }

    def preload_datasets(
        self,
        dataset_page_ids: List[str],
        extract_fn: Callable[..., Dict[str, Set[str]]],
        **extract_kwargs,
    ) -> Dict[str, bool]:
        """
        Preload features for multiple datasets.

        Args:
            dataset_page_ids: List of dataset page IDs to preload
            extract_fn: Function to extract features (e.g., extract_dataset_features_by_type)
            **extract_kwargs: Additional keyword arguments to pass to extract_fn

        Returns:
            Dictionary mapping dataset_page_id → success (bool)
        """
        results: Dict[str, bool] = {}
        logger.info(
            "[CACHE] Preloading features for %d datasets",
            len(dataset_page_ids),
        )

        for dataset_page_id in dataset_page_ids:
            try:
                # Check if already cached and fresh
                cached = self.get_features(dataset_page_id)
                if cached:
                    logger.debug(
                        "[CACHE] Dataset %s already cached, skipping preload",
                        dataset_page_id,
                    )
                    results[dataset_page_id] = True
                    continue

                # Extract and cache
                features_by_type = extract_fn(
                    dataset_page_id=dataset_page_id, **extract_kwargs
                )
                if features_by_type:
                    self.set_features(
                        dataset_page_id=dataset_page_id,
                        features_by_type=features_by_type,
                    )
                    results[dataset_page_id] = True
                else:
                    results[dataset_page_id] = False
            except Exception as e:
                logger.warning(
                    "[CACHE] Error preloading dataset %s: %r",
                    dataset_page_id,
                    e,
                )
                results[dataset_page_id] = False

        success_count = sum(1 for v in results.values() if v)
        logger.info(
            "[CACHE] Preloaded %d/%d datasets successfully",
            success_count,
            len(dataset_page_ids),
        )

        return results


# Global cache instance (singleton pattern)
_feature_cache: Optional[DatasetFeatureCache] = None


def get_feature_cache(ttl_seconds: int = 3600) -> DatasetFeatureCache:
    """
    Get or create the global feature cache instance.

    Args:
        ttl_seconds: Time-to-live for cache entries (only used on first creation)

    Returns:
        Global DatasetFeatureCache instance
    """
    global _feature_cache

    if _feature_cache is None:
        _feature_cache = DatasetFeatureCache(ttl_seconds=ttl_seconds)
        logger.debug(
            "[CACHE] Created global feature cache instance",
        )

    return _feature_cache


def clear_feature_cache(dataset_page_id: Optional[str] = None) -> None:
    """
    Clear the global feature cache.

    Args:
        dataset_page_id: Optional specific dataset to clear. If None, clears entire cache.
    """
    cache = get_feature_cache()
    cache.clear(dataset_page_id=dataset_page_id)

