#!/usr/bin/env python3
"""
Find datasets that have features linked for testing.

Checks datasets to see which ones have features linked to them,
which makes them good candidates for testing signature scoring.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.multi_omics_scoring import (
    extract_dataset_features_by_type,
)
from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
import requests
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def list_all_datasets():
    """List all datasets from Experimental Data Assets database."""
    cfg = get_config()

    if not cfg.notion.exp_data_db_id:
        print("⚠️  Experimental Data Assets database ID not configured")
        return []

    datasets = []
    url = f"{cfg.notion.base_url}/databases/{cfg.notion.exp_data_db_id}/query"

    try:
        has_more = True
        start_cursor = None

        while has_more:
            payload = {"page_size": 100}
            if start_cursor:
                payload["start_cursor"] = start_cursor

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            data = resp.json()
            for page in data.get("results", []):
                props = page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                page_id = page.get("id", "")

                # Get omics type if available
                omics_type_prop = props.get("Omics Type", {}).get("select")
                omics_type = omics_type_prop.get("name", "") if omics_type_prop else None

                datasets.append({
                    "name": name,
                    "page_id": page_id,
                    "omics_type": omics_type,
                })

            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")

        return datasets

    except Exception as e:
        logger.error("Error listing datasets: %r", e)
        return []


def check_dataset_features(dataset_page_id: str, dataset_name: str = "Unknown") -> dict:
    """
    Check how many features are linked to a dataset.

    Returns:
        Dictionary with feature counts by type
    """
    try:
        features = extract_dataset_features_by_type(
            dataset_page_id=dataset_page_id,
            use_cache=False,  # Don't use cache for this check
            force_refresh=True,
        )

        counts = {
            "gene": len(features.get("gene", set())),
            "protein": len(features.get("protein", set())),
            "metabolite": len(features.get("metabolite", set())),
            "lipid": len(features.get("lipid", set())),
        }
        total = sum(counts.values())

        return {
            "page_id": dataset_page_id,
            "name": dataset_name,
            "counts": counts,
            "total": total,
            "has_features": total > 0,
        }

    except Exception as e:
        logger.warning("Error checking features for dataset %s: %r", dataset_page_id, e)
        return {
            "page_id": dataset_page_id,
            "name": dataset_name,
            "counts": {},
            "total": 0,
            "has_features": False,
            "error": str(e),
        }


def main():
    """Main function to find datasets with features."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Find datasets that have features linked"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of datasets to check (default: 20)",
    )
    parser.add_argument(
        "--min-features",
        type=int,
        default=1,
        help="Minimum number of features required (default: 1)",
    )

    args = parser.parse_args()

    print("=" * 70)
    print("🔍 Finding Datasets with Features")
    print("=" * 70)

    # List all datasets
    print("\n📊 Listing all datasets...")
    datasets = list_all_datasets()

    if not datasets:
        print("❌ No datasets found")
        return

    print(f"  Found {len(datasets)} datasets")

    # Check features for each dataset (limit to avoid too many API calls)
    print(f"\n🔍 Checking features for up to {args.limit} datasets...")
    results = []

    for idx, dataset in enumerate(datasets[: args.limit], 1):
        print(f"  [{idx}/{min(len(datasets), args.limit)}] Checking {dataset['name'][:50]}...")
        result = check_dataset_features(
            dataset_page_id=dataset["page_id"],
            dataset_name=dataset["name"],
        )
        results.append(result)

        if result["has_features"]:
            counts = result["counts"]
            print(f"     ✅ {result['total']} features: "
                  f"genes={counts['gene']}, proteins={counts['protein']}, "
                  f"metabolites={counts['metabolite']}, lipids={counts['lipid']}")

    # Filter to datasets with features
    datasets_with_features = [r for r in results if r.get("has_features", False)]

    print("\n" + "=" * 70)
    print("📊 RESULTS")
    print("=" * 70)

    if datasets_with_features:
        print(f"\n✅ Found {len(datasets_with_features)} datasets with features:")
        print()

        for result in datasets_with_features:
            counts = result["counts"]
            print(f"📦 {result['name']}")
            print(f"   ID: {result['page_id']}")
            print(f"   Total features: {result['total']}")
            print(f"   Breakdown: genes={counts['gene']}, proteins={counts['protein']}, "
                  f"metabolites={counts['metabolite']}, lipids={counts['lipid']}")
            print()

        # Show best candidate
        best = max(datasets_with_features, key=lambda x: x["total"])
        print("\n⭐ Best candidate for testing:")
        print(f"   Name: {best['name']}")
        print(f"   ID: {best['page_id']}")
        print(f"   Features: {best['total']} total")

        print("\n💡 Test command:")
        print(f"   python scripts/test_feature_caching.py --dataset-id {best['page_id']} --all-tests")

    else:
        print(f"\n⚠️  No datasets found with features (min: {args.min_features})")
        print(f"   Checked {len(results)} datasets")
        print("\n   Note: This could mean:")
        print("   • Features haven't been linked yet")
        print("   • Need to ingest datasets to create feature links")
        print("   • Feature databases may not be configured")


if __name__ == "__main__":
    main()

