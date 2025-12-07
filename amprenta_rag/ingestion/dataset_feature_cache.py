"""
Dataset feature cache with TTL and LRU eviction.

Reduces repeated external calls by caching per-dataset feature maps.
"""

from __future__ import annotations

import time
from collections import OrderedDict
from threading import Lock
from typing import Any, Callable, Dict, Optional, Set

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

FeatureMap = Dict[str, Set[str]]
CacheEntry = Dict[str, Any]

DEFAULT_TTL_SECONDS = 3600
DEFAULT_MAX_SIZE = 1000


class DatasetFeatureCache:
    """Thread-safe in-memory cache for dataset features."""

    def __init__(
        self,
        ttl_seconds: int = DEFAULT_TTL_SECONDS,
        max_size: int = DEFAULT_MAX_SIZE,
    ):
        self.cache: OrderedDict[str, CacheEntry] = OrderedDict()
        self.default_ttl = ttl_seconds
        self.max_size = max_size
        self._lock = Lock()

        # Stats
        self.hits = 0
        self.misses = 0
        self.evictions = 0

        logger.debug(
            "[FEATURE-CACHE] Initialized (default_ttl=%s, max_size=%s)",
            ttl_seconds,
            max_size,
        )

    def get_features(self, dataset_id: str) -> Optional[FeatureMap]:
        """
        Return cached features for the dataset or None when missing/expired.
        """
        with self._lock:
            entry = self.cache.get(dataset_id)
            if entry is None:
                self.misses += 1
                logger.debug(
                    "[FEATURE-CACHE] Miss for dataset %s (not cached)",
                    dataset_id,
                )
                return None

            if self._is_entry_expired(entry):
                self.cache.pop(dataset_id, None)
                self.evictions += 1
                self.misses += 1
                logger.debug(
                    "[FEATURE-CACHE] Entry expired for dataset %s; evicted",
                    dataset_id,
                )
                return None

            self.cache.move_to_end(dataset_id)
            self.hits += 1
            logger.debug(
                "[FEATURE-CACHE] Hit for dataset %s",
                dataset_id,
            )
            return entry["features"]

    def set_features(
        self,
        dataset_id: str,
        features: FeatureMap,
        ttl: int = DEFAULT_TTL_SECONDS,
    ) -> None:
        """
        Store features for a dataset with optional per-entry TTL override.
        """
        normalized_features: FeatureMap = {k: set(v) for k, v in (features or {}).items()}
        effective_ttl = ttl
        if ttl == DEFAULT_TTL_SECONDS and self.default_ttl != DEFAULT_TTL_SECONDS:
            effective_ttl = self.default_ttl

        timestamp = time.time()

        with self._lock:
            if dataset_id in self.cache:
                self.cache.pop(dataset_id)

            self.cache[dataset_id] = {
                "features": normalized_features,
                "timestamp": timestamp,
                "ttl": effective_ttl,
            }
            self.cache.move_to_end(dataset_id)
            self._evict_if_needed_unlocked()

            logger.debug(
                "[FEATURE-CACHE] Stored dataset %s (types=%d, ttl=%s)",
                dataset_id,
                len(normalized_features),
                effective_ttl,
            )

    def invalidate(self, dataset_id: str) -> None:
        """Remove a dataset from the cache."""
        with self._lock:
            if dataset_id in self.cache:
                self.cache.pop(dataset_id, None)
                self.evictions += 1
                logger.debug(
                    "[FEATURE-CACHE] Invalidated dataset %s",
                    dataset_id,
                )
            else:
                logger.debug(
                    "[FEATURE-CACHE] No entry to invalidate for dataset %s",
                    dataset_id,
                )

    def clear(self) -> None:
        """Clear all cache entries."""
        with self._lock:
            size = len(self.cache)
            self.cache.clear()
            logger.debug(
                "[FEATURE-CACHE] Cleared cache (entries=%s)",
                size,
            )

    def get_stats(self) -> Dict[str, int]:
        """Return cache statistics."""
        with self._lock:
            total_features = sum(
                len(features)
                for entry in self.cache.values()
                for features in entry["features"].values()
            )
            cached_datasets = len(self.cache)

            return {
                "cached_datasets": cached_datasets,
                "size": cached_datasets,
                "total_features": total_features,
                "hits": self.hits,
                "misses": self.misses,
                "evictions": self.evictions,
                "max_size": self.max_size,
            }

    def is_expired(self, dataset_id: str) -> bool:
        """Check if a dataset entry is expired (or missing)."""
        with self._lock:
            entry = self.cache.get(dataset_id)
            if entry is None:
                return True
            return self._is_entry_expired(entry)

    def preload_datasets(
        self,
        dataset_page_ids: list[str],
        extract_fn: Callable[..., FeatureMap],
        **extract_kwargs: Any,
    ) -> Dict[str, bool]:
        """
        Warm the cache by extracting features for the provided dataset IDs.
        """
        results: Dict[str, bool] = {}
        for dataset_id in dataset_page_ids:
            cached = self.get_features(dataset_id)
            if cached is not None:
                results[dataset_id] = True
                continue

            try:
                features = extract_fn(dataset_page_id=dataset_id, **extract_kwargs)
                if features:
                    self.set_features(dataset_id, features)
                    results[dataset_id] = True
                else:
                    results[dataset_id] = False
            except Exception as exc:
                logger.warning(
                    "[FEATURE-CACHE] Error preloading dataset %s: %r",
                    dataset_id,
                    exc,
                )
                results[dataset_id] = False

        return results

    def _evict_if_needed_unlocked(self) -> None:
        if self.max_size <= 0:
            return

        while len(self.cache) > self.max_size:
            evicted_id, _ = self.cache.popitem(last=False)
            self.evictions += 1
            logger.debug(
                "[FEATURE-CACHE] Evicted dataset %s due to max size (%s)",
                evicted_id,
                self.max_size,
            )

    def _is_entry_expired(self, entry: CacheEntry) -> bool:
        ttl_seconds = entry.get("ttl", self.default_ttl)
        timestamp = entry.get("timestamp", 0)
        ts_seconds = (
            timestamp.timestamp()
            if hasattr(timestamp, "timestamp")
            else float(timestamp)
        )
        age = time.time() - ts_seconds
        return age >= ttl_seconds


# Global cache instance (singleton)
_feature_cache: Optional[DatasetFeatureCache] = None
_feature_cache_lock = Lock()


def get_feature_cache() -> DatasetFeatureCache:
    """Return the shared feature cache instance."""
    global _feature_cache

    if _feature_cache is None:
        with _feature_cache_lock:
            if _feature_cache is None:
                _feature_cache = DatasetFeatureCache()
                logger.debug("[FEATURE-CACHE] Created singleton cache instance")

    return _feature_cache


def clear_feature_cache(dataset_id: Optional[str] = None) -> None:
    """Convenience helper to clear or invalidate the shared cache."""
    cache = get_feature_cache()
    if dataset_id:
        cache.invalidate(dataset_id)
    else:
        cache.clear()

