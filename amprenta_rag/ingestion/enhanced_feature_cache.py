"""
Enhanced dataset feature cache with advanced performance optimizations.

This module provides an enhanced caching system with:
- Configurable TTL per dataset type
- LRU eviction policy
- File-based persistence
- Parallel pre-loading
- Performance profiling
"""

from __future__ import annotations

import json
import os
import pickle
from collections import OrderedDict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.performance import timer, get_performance_metrics

logger = get_logger(__name__)


class EnhancedDatasetFeatureCache:
    """
    Enhanced in-memory cache for dataset feature sets with advanced optimizations.
    
    Features:
    - Configurable TTL per dataset type
    - LRU eviction when memory limit reached
    - File-based persistence for cold starts
    - Parallel pre-loading with workers
    - Performance metrics tracking
    """
    
    def __init__(
        self,
        ttl_seconds: int = 3600,
        max_size: Optional[int] = None,
        enable_persistence: bool = True,
        cache_dir: Optional[str] = None,
    ):
        """
        Initialize the enhanced cache.
        
        Args:
            ttl_seconds: Default time-to-live for cache entries (default: 1 hour)
            max_size: Maximum number of datasets to cache (None = unlimited)
            enable_persistence: Enable file-based persistence (default: True)
            cache_dir: Directory for cache files (default: .cache/feature_cache)
        """
        # LRU cache using OrderedDict
        self.cache: OrderedDict[str, Dict[str, Any]] = OrderedDict()
        self.default_ttl = timedelta(seconds=ttl_seconds)
        self.max_size = max_size
        
        # Performance metrics
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        
        # File persistence
        self.enable_persistence = enable_persistence
        if cache_dir:
            self.cache_dir = Path(cache_dir)
        else:
            self.cache_dir = Path(".cache") / "feature_cache"
        
        if self.enable_persistence:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self._load_from_disk()
        
        # TTL configuration per dataset type (can be customized)
        self.ttl_config: Dict[str, timedelta] = {
            "gene": timedelta(seconds=7200),  # 2 hours
            "protein": timedelta(seconds=7200),
            "metabolite": timedelta(seconds=7200),
            "lipid": timedelta(seconds=7200),
        }
        
        logger.info(
            "[CACHE][ENHANCED] Initialized with TTL=%ds, max_size=%s, persistence=%s",
            ttl_seconds,
            max_size or "unlimited",
            enable_persistence,
        )
    
    def _get_ttl(self, omics_type: Optional[str] = None) -> timedelta:
        """Get TTL for a specific omics type or use default."""
        if omics_type and omics_type.lower() in self.ttl_config:
            return self.ttl_config[omics_type.lower()]
        return self.default_ttl
    
    def _evict_lru(self):
        """Evict least recently used entry if cache is full."""
        if self.max_size and len(self.cache) >= self.max_size:
            # Remove oldest entry (first in OrderedDict)
            oldest_key = next(iter(self.cache))
            del self.cache[oldest_key]
            self.evictions += 1
            logger.debug(
                "[CACHE][ENHANCED] Evicted LRU entry: %s",
                oldest_key,
            )
    
    def _move_to_end(self, key: str):
        """Move key to end (most recently used)."""
        if key in self.cache:
            self.cache.move_to_end(key)
    
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
            Dictionary mapping feature_type → set of feature names, or None
        """
        with timer(f"cache_get_{dataset_page_id[:8]}", log_threshold=0.5):
            if force_refresh:
                logger.debug(
                    "[CACHE][ENHANCED] Force refresh requested for dataset %s",
                    dataset_page_id,
                )
                self.misses += 1
                return None
            
            if dataset_page_id not in self.cache:
                logger.debug(
                    "[CACHE][ENHANCED] Cache miss for dataset %s",
                    dataset_page_id,
                )
                self.misses += 1
                return None
            
            entry = self.cache[dataset_page_id]
            age = datetime.now() - entry["timestamp"]
            omics_type = entry.get("omics_type")
            ttl = self._get_ttl(omics_type)
            
            if age > ttl:
                logger.debug(
                    "[CACHE][ENHANCED] Cache entry stale for dataset %s (age: %s, TTL: %s)",
                    dataset_page_id,
                    age,
                    ttl,
                )
                self.misses += 1
                # Remove stale entry
                del self.cache[dataset_page_id]
                return None
            
            # Cache hit - move to end (most recently used)
            self._move_to_end(dataset_page_id)
            self.hits += 1
            
            logger.debug(
                "[CACHE][ENHANCED] Cache hit for dataset %s (age: %s)",
                dataset_page_id,
                age,
            )
            return entry["features_by_type"]
    
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
        with timer(f"cache_set_{dataset_page_id[:8]}", log_threshold=0.5):
            # Evict if needed
            if dataset_page_id not in self.cache:
                self._evict_lru()
            
            feature_type_count = len([k for k, v in features_by_type.items() if v])
            total_features = sum(len(v) for v in features_by_type.values())
            
            self.cache[dataset_page_id] = {
                "features_by_type": features_by_type,
                "directions": directions or {},
                "timestamp": datetime.now(),
                "omics_type": omics_type,
            }
            
            # Move to end (most recently used)
            self._move_to_end(dataset_page_id)
            
            logger.debug(
                "[CACHE][ENHANCED] Cached features for dataset %s (%d types, %d features)",
                dataset_page_id,
                feature_type_count,
                total_features,
            )
            
            # Persist to disk if enabled
            if self.enable_persistence:
                self._save_to_disk(dataset_page_id)
    
    def preload_datasets_parallel(
        self,
        dataset_page_ids: List[str],
        extract_fn: Callable[..., Dict[str, Set[str]]],
        max_workers: int = 5,
        **extract_kwargs,
    ) -> Dict[str, bool]:
        """
        Preload features for multiple datasets in parallel.
        
        Args:
            dataset_page_ids: List of dataset page IDs to preload
            extract_fn: Function to extract features
            max_workers: Maximum number of parallel workers (default: 5)
            **extract_kwargs: Additional keyword arguments to pass to extract_fn
            
        Returns:
            Dictionary mapping dataset_page_id → success (bool)
        """
        results: Dict[str, bool] = {}
        logger.info(
            "[CACHE][ENHANCED] Preloading features for %d datasets (parallel, %d workers)",
            len(dataset_page_ids),
            max_workers,
        )
        
        # Filter out already-cached fresh entries
        to_preload = []
        for dataset_page_id in dataset_page_ids:
            cached = self.get_features(dataset_page_id)
            if cached:
                results[dataset_page_id] = True
            else:
                to_preload.append(dataset_page_id)
        
        if not to_preload:
            logger.info(
                "[CACHE][ENHANCED] All %d datasets already cached",
                len(dataset_page_ids),
            )
            return results
        
        logger.info(
            "[CACHE][ENHANCED] Preloading %d datasets (skipped %d already cached)",
            len(to_preload),
            len(dataset_page_ids) - len(to_preload),
        )
        
        def preload_one(dataset_page_id: str) -> tuple[str, bool]:
            """Preload a single dataset."""
            try:
                features_by_type = extract_fn(
                    dataset_page_id=dataset_page_id,
                    use_cache=False,  # Don't use cache during preloading
                    **extract_kwargs,
                )
                if features_by_type:
                    # Determine omics type from features
                    omics_type = None
                    for ft, features in features_by_type.items():
                        if features:
                            omics_type = ft
                            break
                    
                    self.set_features(
                        dataset_page_id=dataset_page_id,
                        features_by_type=features_by_type,
                        omics_type=omics_type,
                    )
                    return (dataset_page_id, True)
                else:
                    return (dataset_page_id, False)
            except Exception as e:
                logger.warning(
                    "[CACHE][ENHANCED] Error preloading dataset %s: %r",
                    dataset_page_id,
                    e,
                )
                return (dataset_page_id, False)
        
        # Parallel pre-loading
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(preload_one, dataset_id): dataset_id
                for dataset_id in to_preload
            }
            
            for future in as_completed(futures):
                dataset_id, success = future.result()
                results[dataset_id] = success
        
        success_count = sum(1 for v in results.values() if v)
        logger.info(
            "[CACHE][ENHANCED] Preloaded %d/%d datasets successfully",
            success_count,
            len(dataset_page_ids),
        )
        
        return results
    
    def _save_to_disk(self, dataset_page_id: str):
        """Save a single cache entry to disk."""
        if dataset_page_id not in self.cache:
            return
        
        try:
            entry = self.cache[dataset_page_id]
            
            # Convert sets to lists for JSON serialization
            cache_data = {
                "features_by_type": {
                    k: list(v) for k, v in entry["features_by_type"].items()
                },
                "directions": entry["directions"],
                "timestamp": entry["timestamp"].isoformat(),
                "omics_type": entry.get("omics_type"),
            }
            
            cache_file = self.cache_dir / f"{dataset_page_id.replace('-', '_')}.json"
            with open(cache_file, "w") as f:
                json.dump(cache_data, f, indent=2)
            
        except Exception as e:
            logger.debug(
                "[CACHE][ENHANCED] Error saving cache to disk for %s: %r",
                dataset_page_id,
                e,
            )
    
    def _load_from_disk(self):
        """Load cache entries from disk."""
        if not self.cache_dir.exists():
            return
        
        logger.info(
            "[CACHE][ENHANCED] Loading cache from disk: %s",
            self.cache_dir,
        )
        
        loaded = 0
        for cache_file in self.cache_dir.glob("*.json"):
            try:
                dataset_page_id = cache_file.stem.replace("_", "-")
                
                with open(cache_file) as f:
                    cache_data = json.load(f)
                
                # Convert lists back to sets
                features_by_type = {
                    k: set(v) for k, v in cache_data["features_by_type"].items()
                }
                
                timestamp = datetime.fromisoformat(cache_data["timestamp"])
                
                # Only load if not stale
                age = datetime.now() - timestamp
                if age <= self.default_ttl:
                    self.cache[dataset_page_id] = {
                        "features_by_type": features_by_type,
                        "directions": cache_data.get("directions", {}),
                        "timestamp": timestamp,
                        "omics_type": cache_data.get("omics_type"),
                    }
                    loaded += 1
                else:
                    # Remove stale file
                    cache_file.unlink()
            
            except Exception as e:
                logger.debug(
                    "[CACHE][ENHANCED] Error loading cache file %s: %r",
                    cache_file,
                    e,
                )
        
        logger.info(
            "[CACHE][ENHANCED] Loaded %d cache entries from disk",
            loaded,
        )
    
    def get_stats(self) -> Dict[str, Any]:
        """Get comprehensive cache statistics."""
        total_features = sum(
            len(features)
            for entry in self.cache.values()
            for features in entry["features_by_type"].values()
        )
        
        hit_rate = (
            self.hits / (self.hits + self.misses) if (self.hits + self.misses) > 0 else 0.0
        )
        
        return {
            "cached_datasets": len(self.cache),
            "total_features": total_features,
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate": f"{hit_rate:.1%}",
            "evictions": self.evictions,
            "max_size": self.max_size,
            "persistence_enabled": self.enable_persistence,
        }
    
    def clear(self, dataset_page_id: Optional[str] = None) -> None:
        """Clear cache entry(ies)."""
        if dataset_page_id:
            if dataset_page_id in self.cache:
                del self.cache[dataset_page_id]
                
                # Remove from disk if persistence enabled
                if self.enable_persistence:
                    cache_file = self.cache_dir / f"{dataset_page_id.replace('-', '_')}.json"
                    if cache_file.exists():
                        cache_file.unlink()
                
                logger.debug(
                    "[CACHE][ENHANCED] Cleared cache for dataset %s",
                    dataset_page_id,
                )
        else:
            count = len(self.cache)
            self.cache.clear()
            
            # Clear disk cache
            if self.enable_persistence and self.cache_dir.exists():
                for cache_file in self.cache_dir.glob("*.json"):
                    cache_file.unlink()
            
            logger.debug(
                "[CACHE][ENHANCED] Cleared entire cache (%d entries)",
                count,
            )


# Global enhanced cache instance
_enhanced_cache: Optional[EnhancedDatasetFeatureCache] = None


def get_enhanced_feature_cache(
    ttl_seconds: int = 3600,
    max_size: Optional[int] = None,
    enable_persistence: bool = True,
) -> EnhancedDatasetFeatureCache:
    """
    Get or create the global enhanced feature cache instance.
    
    Args:
        ttl_seconds: Time-to-live for cache entries (only used on first creation)
        max_size: Maximum cache size (only used on first creation)
        enable_persistence: Enable file persistence (only used on first creation)
        
    Returns:
        Global EnhancedDatasetFeatureCache instance
    """
    global _enhanced_cache
    
    if _enhanced_cache is None:
        _enhanced_cache = EnhancedDatasetFeatureCache(
            ttl_seconds=ttl_seconds,
            max_size=max_size,
            enable_persistence=enable_persistence,
        )
    
    return _enhanced_cache

