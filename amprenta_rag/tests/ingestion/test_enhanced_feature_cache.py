from __future__ import annotations

import json
from pathlib import Path

from amprenta_rag.ingestion.enhanced_feature_cache import EnhancedDatasetFeatureCache


def test_set_and_get_features_with_ttl(monkeypatch):
    cache = EnhancedDatasetFeatureCache(ttl_seconds=1, max_size=2, enable_persistence=False)

    cache.set_features("ds1", {"gene": {"A", "B"}})
    assert cache.get_features("ds1") == {"gene": {"A", "B"}}
    assert cache.hits == 1
    assert cache.misses == 0

    cache.cache["ds1"]["timestamp"] = cache.cache["ds1"]["timestamp"] - cache.default_ttl - cache.default_ttl
    assert cache.get_features("ds1") is None
    assert cache.misses >= 1


def test_persistence_round_trip(tmp_path: Path):
    cache_dir = tmp_path / "cache"
    cache = EnhancedDatasetFeatureCache(ttl_seconds=3600, enable_persistence=True, cache_dir=str(cache_dir))

    cache.set_features("ds1", {"gene": {"A"}}, directions={"gene": {"A": "up"}}, omics_type="gene")

    files = list(cache_dir.glob("*.json"))
    assert files
    data = json.loads(files[0].read_text())
    assert data["features_by_type"]["gene"] == ["A"]

    cache2 = EnhancedDatasetFeatureCache(ttl_seconds=3600, enable_persistence=True, cache_dir=str(cache_dir))
    assert cache2.get_features("ds1") == {"gene": {"A"}}


def test_preload_datasets_parallel_skips_cached(monkeypatch):
    cache = EnhancedDatasetFeatureCache(enable_persistence=False)

    cache.set_features("ds1", {"gene": {"A"}})

    def fake_extract(dataset_page_id, use_cache=False):
        return {"protein": {"P1"}} if dataset_page_id == "ds2" else {}

    results = cache.preload_datasets_parallel(["ds1", "ds2"], extract_fn=fake_extract, max_workers=2)
    assert results["ds1"] is True
    assert results["ds2"] is True
    assert cache.get_features("ds2") == {"protein": {"P1"}}
from __future__ import annotations

import json
from pathlib import Path

from amprenta_rag.ingestion.enhanced_feature_cache import EnhancedDatasetFeatureCache


def test_set_and_get_features_with_ttl(monkeypatch):
    cache = EnhancedDatasetFeatureCache(ttl_seconds=1, max_size=2, enable_persistence=False)

    cache.set_features("ds1", {"gene": {"A", "B"}})
    assert cache.get_features("ds1") == {"gene": {"A", "B"}}
    assert cache.hits == 1
    assert cache.misses == 0

    # Mark entry stale by manipulating timestamp
    cache.cache["ds1"]["timestamp"] = cache.cache["ds1"]["timestamp"] - cache.default_ttl - cache.default_ttl
    assert cache.get_features("ds1") is None
    assert cache.misses >= 1


def test_persistence_round_trip(tmp_path: Path):
    cache_dir = tmp_path / "cache"
    cache = EnhancedDatasetFeatureCache(ttl_seconds=3600, enable_persistence=True, cache_dir=str(cache_dir))

    cache.set_features("ds1", {"gene": {"A"}}, directions={"gene": {"A": "up"}}, omics_type="gene")

    files = list(cache_dir.glob("*.json"))
    assert files
    data = json.loads(files[0].read_text())
    assert data["features_by_type"]["gene"] == ["A"]

    cache2 = EnhancedDatasetFeatureCache(ttl_seconds=3600, enable_persistence=True, cache_dir=str(cache_dir))
    assert cache2.get_features("ds1") == {"gene": {"A"}}


def test_preload_datasets_parallel_skips_cached(monkeypatch):
    cache = EnhancedDatasetFeatureCache(enable_persistence=False)

    cache.set_features("ds1", {"gene": {"A"}})

    def fake_extract(dataset_page_id, use_cache=False):
        return {"protein": {"P1"}} if dataset_page_id == "ds2" else {}

    results = cache.preload_datasets_parallel(["ds1", "ds2"], extract_fn=fake_extract, max_workers=2)
    assert results["ds1"] is True
    assert results["ds2"] is True
    assert cache.get_features("ds2") == {"protein": {"P1"}}

