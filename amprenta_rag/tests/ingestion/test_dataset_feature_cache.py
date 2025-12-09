"""
Tests for dataset feature cache.
"""

import threading
import time
from typing import List

import pytest

from amprenta_rag.ingestion.dataset_feature_cache import (
    DatasetFeatureCache,
    clear_feature_cache,
    get_feature_cache,
)


def _features() -> dict[str, set[str]]:
    return {
        "gene": {"TP53", "TNF"},
        "protein": {"P04637"},
        "metabolite": set(),
        "lipid": set(),
    }


def test_cache_set_and_get():
    cache = DatasetFeatureCache()
    features = _features()

    cache.set_features("test-dataset-id", features)
    retrieved = cache.get_features("test-dataset-id")

    assert retrieved == features
    stats = cache.get_stats()
    assert stats["hits"] == 1
    assert stats["misses"] == 0
    assert stats["cached_datasets"] == 1
    assert stats["size"] == 1
    assert stats["total_features"] == sum(len(v) for v in features.values())


def test_cache_miss_and_expiry():
    cache = DatasetFeatureCache(ttl_seconds=1)
    features = _features()

    # Initial miss
    assert cache.get_features("missing-id") is None
    stats = cache.get_stats()
    assert stats["hits"] == 0
    assert stats["misses"] == 1

    # Expired entry behaves like miss and is removed
    cache.set_features("test-id", features)
    cache.cache["test-id"]["timestamp"] = time.time() - 2

    assert cache.get_features("test-id") is None
    assert cache.is_expired("test-id")
    stats = cache.get_stats()
    # 2 total misses: one for missing-id, one for expired test-id
    assert stats["misses"] == 2
    # Expired entry triggers an eviction
    assert stats["evictions"] >= 1


def test_invalidate_and_clear():
    cache = DatasetFeatureCache()
    cache.set_features("id1", _features())
    cache.set_features("id2", _features())

    cache.invalidate("id1")
    assert "id1" not in cache.cache

    # Invalidate should increment evictions
    stats = cache.get_stats()
    assert stats["evictions"] >= 1

    cache.clear()
    assert len(cache.cache) == 0
    stats = cache.get_stats()
    assert stats["cached_datasets"] == 0
    assert stats["size"] == 0
    assert stats["total_features"] == 0


def test_lru_eviction_tracks_stats():
    cache = DatasetFeatureCache(max_size=2)
    cache.set_features("id1", _features())
    cache.set_features("id2", _features())
    cache.set_features("id3", _features())

    # Oldest entry (id1) should be evicted
    assert "id1" not in cache.cache
    stats = cache.get_stats()
    assert stats["evictions"] >= 1
    assert stats["cached_datasets"] == 2
    assert stats["size"] == 2


def test_cache_hit_and_miss_stats():
    cache = DatasetFeatureCache()
    cache.set_features("id1", _features())

    # One hit
    assert cache.get_features("id1") is not None
    # One miss
    assert cache.get_features("missing-id") is None

    stats = cache.get_stats()
    assert stats["hits"] == 1
    assert stats["misses"] == 1


def test_singleton_cache():
    clear_feature_cache()
    cache1 = get_feature_cache()
    cache2 = get_feature_cache()

    assert cache1 is cache2


@pytest.mark.parametrize("ttl_seconds,expected_expired", [(0, True), (1, False)])
def test_is_expired_edge_cases(ttl_seconds: int, expected_expired: bool):
    cache = DatasetFeatureCache(ttl_seconds=ttl_seconds)
    cache.set_features("id", _features())

    # Manually adjust timestamp to now to control expiry behavior
    cache.cache["id"]["timestamp"] = time.time()

    assert cache.is_expired("id") is expected_expired


def test_thread_safety_under_concurrent_access():
    """
    Concurrent get/set operations should not corrupt the cache or stats.
    """
    cache = DatasetFeatureCache(max_size=10)
    dataset_ids = [f"id{i}" for i in range(5)]
    features = _features()

    def writer():
        for _ in range(100):
            for did in dataset_ids:
                cache.set_features(did, features)

    def reader():
        for _ in range(100):
            for did in dataset_ids:
                _ = cache.get_features(did)

    threads: List[threading.Thread] = []
    for _ in range(5):
        threads.append(threading.Thread(target=writer))
        threads.append(threading.Thread(target=reader))

    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # Cache should not exceed max_size and should not raise exceptions
    stats = cache.get_stats()
    assert stats["cached_datasets"] <= cache.max_size
    # Hits + misses should be non-negative and consistent
    assert stats["hits"] >= 0
    assert stats["misses"] >= 0

