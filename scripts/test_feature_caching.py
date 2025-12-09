#!/usr/bin/env python3
"""
Test script for feature caching implementation.

Tests the caching system with real datasets to verify:
1. Cache hit/miss behavior
2. Performance improvements
3. Batch scoring functionality
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.dataset_feature_cache import (
    clear_feature_cache,
)
from amprenta_rag.ingestion.multi_omics_scoring import (
    extract_dataset_features_by_type,
)
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset,
)
from amprenta_rag.ingestion.batch_signature_scoring import (
    score_datasets_against_signatures_batch,
    get_cache_stats,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def test_cache_hit_miss(dataset_page_id: str) -> None:
    """
    Test cache hit/miss behavior for a single dataset.

    Args:
        dataset_page_id: Notion page ID of dataset to test
    """
    print("\n" + "=" * 70)
    print("🧪 Test 1: Cache Hit/Miss Behavior")
    print("=" * 70)

    # Clear cache first
    clear_feature_cache(dataset_page_id)

    # First call - should miss cache and hit Notion
    print(f"\n📊 First call (cache miss expected)...")
    start_time = time.time()
    features1 = extract_dataset_features_by_type(
        dataset_page_id=dataset_page_id,
        use_cache=True,
        force_refresh=False,
    )
    time1 = time.time() - start_time
    print(f"  ⏱️  Time: {time1:.3f}s")
    print(f"  📦 Features extracted: {sum(len(v) for v in features1.values())}")

    # Second call - should hit cache
    print(f"\n📊 Second call (cache hit expected)...")
    start_time = time.time()
    features2 = extract_dataset_features_by_type(
        dataset_page_id=dataset_page_id,
        use_cache=True,
        force_refresh=False,
    )
    time2 = time.time() - start_time
    print(f"  ⏱️  Time: {time2:.3f}s")
    print(f"  📦 Features extracted: {sum(len(v) for v in features2.values())}")

    # Verify results match
    assert features1 == features2, "Features should match between calls"
    print(f"  ✅ Features match between calls")

    # Calculate speedup
    if time1 > 0 and time2 > 0:
        speedup = time1 / time2
        print(f"\n⚡ Performance Improvement:")
        print(f"  First call (cache miss): {time1:.3f}s")
        print(f"  Second call (cache hit):  {time2:.3f}s")
        print(f"  Speedup: {speedup:.1f}x faster with cache!")

    # Check cache stats
    stats = get_cache_stats()
    print(f"\n📊 Cache Statistics:")
    print(f"  Cached datasets: {stats.get('cached_datasets', 0)}")
    print(f"  Total features: {stats.get('total_features', 0)}")


def test_signature_scoring_with_cache(dataset_page_id: str) -> None:
    """
    Test signature scoring with caching enabled.

    Args:
        dataset_page_id: Notion page ID of dataset to test
    """
    print("\n" + "=" * 70)
    print("🧪 Test 2: Signature Scoring with Cache")
    print("=" * 70)

    # Clear cache first
    clear_feature_cache(dataset_page_id)

    # Score with cache (first call - will populate cache)
    print(f"\n📊 First scoring (will populate cache)...")
    start_time = time.time()
    matches1 = find_matching_signatures_for_dataset(
        dataset_page_id=dataset_page_id,
        overlap_threshold=0.3,
    )
    time1 = time.time() - start_time
    print(f"  ⏱️  Time: {time1:.3f}s")
    print(f"  🎯 Matches found: {len(matches1)}")

    # Score again with cache (should use cache)
    print(f"\n📊 Second scoring (should use cache)...")
    start_time = time.time()
    matches2 = find_matching_signatures_for_dataset(
        dataset_page_id=dataset_page_id,
        overlap_threshold=0.3,
    )
    time2 = time.time() - start_time
    print(f"  ⏱️  Time: {time2:.3f}s")
    print(f"  🎯 Matches found: {len(matches2)}")

    # Verify results match
    assert len(matches1) == len(matches2), "Match counts should be the same"
    print(f"  ✅ Results match between calls")

    # Calculate speedup
    if time1 > 0 and time2 > 0:
        speedup = time1 / time2
        print(f"\n⚡ Performance Improvement:")
        print(f"  First scoring: {time1:.3f}s")
        print(f"  Second scoring: {time2:.3f}s")
        print(f"  Speedup: {speedup:.1f}x faster with cache!")


def test_batch_scoring(dataset_page_ids: list[str]) -> None:
    """
    Test batch scoring with cache pre-loading.

    Args:
        dataset_page_ids: List of dataset page IDs to test
    """
    print("\n" + "=" * 70)
    print("🧪 Test 3: Batch Scoring with Cache Pre-loading")
    print("=" * 70)

    # Clear cache
    clear_feature_cache()

    print(f"\n📊 Testing batch scoring for {len(dataset_page_ids)} datasets...")

    # Batch score with cache pre-loading
    start_time = time.time()
    results = score_datasets_against_signatures_batch(
        dataset_page_ids=dataset_page_ids,
        overlap_threshold=0.3,
        preload_cache=True,
        use_cache=True,
    )
    total_time = time.time() - start_time

    print(f"\n⏱️  Total time: {total_time:.3f}s")
    print(f"  Average per dataset: {total_time/len(dataset_page_ids):.3f}s")

    # Show results
    print(f"\n📊 Results:")
    for dataset_id, matches in results.items():
        print(f"  Dataset {dataset_id[:8]}...: {len(matches)} matches")

    # Cache stats
    stats = get_cache_stats()
    print(f"\n📊 Cache Statistics:")
    print(f"  Cached datasets: {stats.get('cached_datasets', 0)}")
    print(f"  Total features: {stats.get('total_features', 0)}")


def main() -> None:
    """Main test function."""
    import argparse

    parser = argparse.ArgumentParser(description="Test feature caching implementation")
    parser.add_argument(
        "--dataset-id",
        type=str,
        help="Single dataset page ID to test (with dashes)",
    )
    parser.add_argument(
        "--batch-datasets",
        nargs="+",
        help="Multiple dataset page IDs for batch testing",
    )
    parser.add_argument(
        "--all-tests",
        action="store_true",
        help="Run all tests if dataset-id provided",
    )
    parser.add_argument(
        "--clear-cache",
        action="store_true",
        help="Clear cache before testing",
    )

    args = parser.parse_args()

    print("=" * 70)
    print("🧪 Feature Caching Test Suite")
    print("=" * 70)

    # Clear cache if requested
    if args.clear_cache:
        print("\n🗑️  Clearing cache...")
        clear_feature_cache()
        print("  ✅ Cache cleared")

    # Single dataset tests
    if args.dataset_id:
        if args.all_tests:
            test_cache_hit_miss(args.dataset_id)
            test_signature_scoring_with_cache(args.dataset_id)
        else:
            # Default: just test cache hit/miss
            test_cache_hit_miss(args.dataset_id)

    # Batch tests
    if args.batch_datasets:
        test_batch_scoring(args.batch_datasets)

    if not args.dataset_id and not args.batch_datasets:
        print("\n❌ No dataset IDs provided!")
        print("\nUsage examples:")
        print("  python scripts/test_feature_caching.py --dataset-id <page-id>")
        print(
            "  python scripts/test_feature_caching.py --dataset-id <page-id> --all-tests"
        )
        print(
            "  python scripts/test_feature_caching.py --batch-datasets <id1> <id2> <id3>"
        )
        sys.exit(1)

    print("\n" + "=" * 70)
    print("✅ Testing Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()

