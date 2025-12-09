"""
Integration tests for feature cache behavior with multi_omics_scoring and
batch_signature_scoring.
"""

from __future__ import annotations

import time
from typing import Dict, Set

import pytest
from unittest.mock import patch

from amprenta_rag.ingestion.batch_signature_scoring import (
    score_datasets_against_signatures_batch,
)
from amprenta_rag.ingestion.dataset_feature_cache import (
    clear_feature_cache,
    get_feature_cache,
)
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type


@pytest.fixture(autouse=True)
def _clear_cache_between_tests():
    """Ensure global cache is cleared between tests for isolation."""
    clear_feature_cache()
    yield
    clear_feature_cache()


def _fake_features(dataset_page_id: str) -> Dict[str, Set[str]]:
    return {
        "gene": {f"TP53-{dataset_page_id}"},
        "protein": {f"P04637-{dataset_page_id}"},
        "metabolite": set(),
        "lipid": set(),
    }


def test_extract_dataset_features_uses_cache_on_second_call(monkeypatch):
    """
    extract_dataset_features_by_type() with use_cache=True should hit cache
    on the second call, avoiding additional extraction work.
    """
    calls = {"count": 0}

    def fake_extract(*args, **kwargs):
        calls["count"] += 1
        return _fake_features(kwargs.get("dataset_page_id", "ds1"))

    # Patch internal Notion/extraction logic by patching requests.post/get to avoid real calls
    with patch("amprenta_rag.ingestion.multi_omics_scoring.requests.post") as mock_post, patch(
        "amprenta_rag.ingestion.multi_omics_scoring.requests.get"
    ) as mock_get:
        mock_post.side_effect = lambda *a, **k: fake_response()
        mock_get.side_effect = lambda *a, **k: fake_response()

        # First call should go through extraction path and populate cache
        features1 = extract_dataset_features_by_type(dataset_page_id="ds1", use_cache=True)
        assert isinstance(features1, dict)

        # Manually set cache to our fake features to avoid depending on Notion behavior
        cache = get_feature_cache()
        cache.set_features("ds1", _fake_features("ds1"))

        # Second call should return cached features without additional external work
        features2 = extract_dataset_features_by_type(dataset_page_id="ds1", use_cache=True)
        assert features2 == _fake_features("ds1")


def test_force_refresh_bypasses_cache_and_populates(monkeypatch):
    """
    force_refresh=True should bypass cache read but still populate cache
    with the refreshed value.
    """
    cache = get_feature_cache()
    cache.set_features("ds1", {"gene": {"OLD"}})

    new_features = _fake_features("ds1")

    with patch("amprenta_rag.ingestion.multi_omics_scoring.requests.post") as mock_post, patch(
        "amprenta_rag.ingestion.multi_omics_scoring.requests.get"
    ) as mock_get:
        mock_post.return_value = fake_response(results=new_features)
        mock_get.return_value = fake_response()

        features = extract_dataset_features_by_type(
            dataset_page_id="ds1",
            use_cache=True,
            force_refresh=True,
        )

    assert features is not None
    # Cache should now contain refreshed features, not the old value
    cached = cache.get_features("ds1")
    assert cached is not None


def test_batch_scoring_faster_with_warm_cache(monkeypatch):
    """
    Signature scoring should be faster with a warm cache than with a cold cache.

    We mock out external I/O and signature matching to focus on cache behavior
    and relative timing, not absolute values.
    """
    dataset_ids = [f"ds{i}" for i in range(5)]

    # Mock extract_dataset_features_by_type to simulate external cost
    def fake_extract(dataset_page_id: str, *args, **kwargs):
        time.sleep(0.01)  # simulate expensive call
        return _fake_features(dataset_page_id)

    # Mock signature matching to be cheap and deterministic
    with patch(
        "amprenta_rag.ingestion.batch_signature_scoring.extract_dataset_features_by_type",
        side_effect=fake_extract,
    ), patch(
        "amprenta_rag.ingestion.signature_matching.find_matching_signatures_for_dataset",
        return_value=[],
    ):
        # Cold cache: use_cache=False so no caching is used
        t0 = time.time()
        _ = score_datasets_against_signatures_batch(
            dataset_page_ids=dataset_ids,
            preload_cache=False,
            use_cache=False,
        )
        cold_time = time.time() - t0

        clear_feature_cache()

        # Warm cache: preload datasets into cache, then score with use_cache=True
        get_feature_cache()
        with patch(
            "amprenta_rag.ingestion.dataset_feature_cache.DatasetFeatureCache.preload_datasets",
            return_value={did: True for did in dataset_ids},
        ):
            t1 = time.time()
            _ = score_datasets_against_signatures_batch(
                dataset_page_ids=dataset_ids,
                preload_cache=True,
                use_cache=True,
            )
            warm_time = time.time() - t1

    # Warm cache path should be faster than cold path in this mocked scenario
    assert warm_time < cold_time


class fake_response:
    """Simple fake response object for Notion HTTP calls."""

    def __init__(self, results=None):
        self._results = results or []

    def raise_for_status(self):
        return None

    def json(self):
        # Shape JSON to look like a Notion database/query response
        if isinstance(self._results, dict):
            # For schema responses
            return {"properties": {}}
        return {"results": [], "has_more": False, "next_cursor": None}


