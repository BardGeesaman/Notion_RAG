from __future__ import annotations

from datetime import datetime, timedelta

from amprenta_rag.ingestion.enhanced_feature_cache import (
    EnhancedDatasetFeatureCache,
    get_enhanced_feature_cache,
)


def _features(dataset_id: str) -> dict[str, set[str]]:
    return {"gene": {f"{dataset_id}-g1"}, "protein": {f"{dataset_id}-p1"}}


def test_set_and_get_features_hits_and_misses():
    cache = EnhancedDatasetFeatureCache(enable_persistence=False)

    # Miss before set
    assert cache.get_features("d1") is None
    assert cache.misses == 1

    cache.set_features("d1", _features("d1"), omics_type="gene")
    res = cache.get_features("d1")
    assert res["gene"] == {"d1-g1"}
    assert cache.hits == 1

    # Force refresh triggers miss
    assert cache.get_features("d1", force_refresh=True) is None
    assert cache.misses == 2


def test_ttl_expiry_marks_entry_stale():
    cache = EnhancedDatasetFeatureCache(ttl_seconds=1, enable_persistence=False)
    cache.set_features("d1", _features("d1"), omics_type="gene")

    # Age the entry beyond TTL
    cache.cache["d1"]["timestamp"] = datetime.now() - timedelta(seconds=5)
    assert cache.get_features("d1") is None
    assert cache.misses >= 1
    assert "d1" not in cache.cache


def test_lru_eviction_when_max_size_exceeded():
    cache = EnhancedDatasetFeatureCache(max_size=1, enable_persistence=False)
    cache.set_features("d1", _features("d1"))
    cache.set_features("d2", _features("d2"))

    assert "d1" not in cache.cache  # evicted
    assert cache.evictions == 1
    assert cache.get_features("d2") is not None


def test_preload_datasets_parallel_respects_cached(monkeypatch):
    cache = EnhancedDatasetFeatureCache(enable_persistence=False)

    # Preload should skip already cached
    cache.set_features("d1", _features("d1"))

    def fake_extract(dataset_page_id: str, use_cache: bool = False, **kwargs):
        if dataset_page_id == "d2":
            return _features("d2")
        return {}

    results = cache.preload_datasets_parallel(["d1", "d2", "d3"], fake_extract, max_workers=2)
    assert results["d1"] is True  # already cached
    assert results["d2"] is True  # loaded via extract
    assert results["d3"] is False  # extract returned empty


def test_persistence_round_trip(tmp_path):
    cache_dir = tmp_path / "cache"
    cache = EnhancedDatasetFeatureCache(enable_persistence=True, cache_dir=str(cache_dir))
    cache.set_features("d1", _features("d1"), omics_type="gene")

    # New cache instance should load from disk
    cache2 = EnhancedDatasetFeatureCache(enable_persistence=True, cache_dir=str(cache_dir))
    loaded = cache2.get_features("d1")
    assert loaded is not None
    assert loaded["gene"] == {"d1-g1"}


def test_get_stats_and_clear():
    cache = EnhancedDatasetFeatureCache(enable_persistence=False)
    cache.set_features("d1", _features("d1"))
    cache.set_features("d2", _features("d2"))

    stats = cache.get_stats()
    assert stats["cached_datasets"] == 2
    assert stats["total_features"] >= 2

    cache.clear("d1")
    assert cache.get_features("d1") is None

    cache.clear()
    assert cache.cache == {}


def test_global_singleton_reuse():
    c1 = get_enhanced_feature_cache(enable_persistence=False)
    c2 = get_enhanced_feature_cache(enable_persistence=False)
    assert c1 is c2

