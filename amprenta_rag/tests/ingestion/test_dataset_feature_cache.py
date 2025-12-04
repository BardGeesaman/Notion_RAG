"""
Tests for dataset feature cache.
"""

import pytest
from datetime import datetime, timedelta

from amprenta_rag.ingestion.dataset_feature_cache import (
    DatasetFeatureCache,
    clear_feature_cache,
    get_feature_cache,
)


def test_cache_set_and_get():
    """Test setting and getting features from cache."""
    cache = DatasetFeatureCache(ttl_seconds=3600)
    features = {
        "gene": {"TP53", "TNF"},
        "protein": {"P04637"},
        "metabolite": set(),
        "lipid": set(),
    }
    cache.set_features("test-dataset-id", features)
    retrieved = cache.get_features("test-dataset-id")
    assert retrieved == features


def test_cache_miss():
    """Test cache miss for non-existent dataset."""
    cache = DatasetFeatureCache()
    result = cache.get_features("non-existent-id")
    assert result is None


def test_cache_ttl_expiry():
    """Test that cache entries expire after TTL."""
    cache = DatasetFeatureCache(ttl_seconds=1)  # 1 second TTL
    features = {"gene": {"TP53"}, "protein": set(), "metabolite": set(), "lipid": set()}
    cache.set_features("test-id", features)
    
    # Should be cached immediately
    assert cache.get_features("test-id") is not None
    
    # Wait for expiry (simulate by manipulating timestamp)
    cache.cache["test-id"]["timestamp"] = datetime.now() - timedelta(seconds=2)
    
    # Should be expired now
    result = cache.get_features("test-id")
    assert result is None


def test_cache_clear():
    """Test clearing cache entries."""
    cache = DatasetFeatureCache()
    cache.set_features("id1", {"gene": {"TP53"}, "protein": set(), "metabolite": set(), "lipid": set()})
    cache.set_features("id2", {"gene": {"TNF"}, "protein": set(), "metabolite": set(), "lipid": set()})
    
    assert len(cache.cache) == 2
    cache.clear("id1")
    assert len(cache.cache) == 1
    assert "id2" in cache.cache
    
    cache.clear()
    assert len(cache.cache) == 0


def test_cache_stats():
    """Test cache statistics."""
    cache = DatasetFeatureCache()
    cache.set_features("id1", {
        "gene": {"TP53", "TNF"},
        "protein": {"P04637"},
        "metabolite": set(),
        "lipid": set(),
    })
    stats = cache.get_stats()
    assert stats["cached_datasets"] == 1
    assert stats["total_features"] == 3  # 2 genes + 1 protein


def test_singleton_cache():
    """Test that get_feature_cache returns singleton."""
    cache1 = get_feature_cache()
    cache2 = get_feature_cache()
    assert cache1 is cache2

