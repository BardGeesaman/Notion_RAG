#!/usr/bin/env python3
"""
Warm the dataset feature cache for faster signature scoring.

Supports warming by:
- Program ID (Notion relation)
- All datasets (Experimental Data Assets DB)
- Explicit dataset IDs
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import requests
from tqdm import tqdm

# Ensure repository imports work when run as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def notion_headers() -> Dict[str, str]:
    """Build Notion headers; fall back to stub if API key missing."""
    cfg = get_config()
    api_key = getattr(cfg.notion, "api_key", "")
    version = getattr(cfg.notion, "version", "2022-06-28")
    base_headers = {"Notion-Version": version, "Content-Type": "application/json"}

    if api_key:
        base_headers["Authorization"] = f"Bearer {api_key}"
    else:
        base_headers["Authorization"] = "Bearer NOTION_REMOVED"

    return base_headers


def _query_datasets(filter_obj: Optional[Dict] = None) -> List[str]:
    """Query Experimental Data Assets DB and return dataset page IDs."""
    cfg = get_config()
    exp_data_db_id = getattr(cfg.notion, "exp_data_db_id", "")
    if not exp_data_db_id:
        logger.warning("[CACHE-WARM] Experimental Data Assets DB ID not configured")
        return []

    url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
    dataset_ids: List[str] = []
    seen: Set[str] = set()

    has_more = True
    start_cursor = None
    while has_more:
        payload: Dict[str, object] = {"page_size": 100}
        if start_cursor:
            payload["start_cursor"] = start_cursor
        if filter_obj:
            payload["filter"] = filter_obj

        try:
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning(
                "[CACHE-WARM] Error querying datasets (cursor=%s): %r",
                start_cursor,
                exc,
            )
            break

        for page in data.get("results", []) or []:
            page_id = page.get("id")
            if page_id and page_id not in seen:
                dataset_ids.append(page_id)
                seen.add(page_id)

        has_more = data.get("has_more", False)
        start_cursor = data.get("next_cursor")

    return dataset_ids


def fetch_all_dataset_ids() -> List[str]:
    """List all dataset page IDs from Experimental Data Assets DB."""
    datasets = _query_datasets()
    if not datasets:
        logger.warning("[CACHE-WARM] No datasets found in Experimental Data Assets DB")
    return datasets


def fetch_dataset_ids_for_program(program_id: str) -> List[str]:
    """List dataset IDs linked to a Program via relation properties."""
    property_candidates = [
        "Related Programs",
        "Programs",
        "Program",
    ]

    collected: List[str] = []
    seen: Set[str] = set()

    for prop in property_candidates:
        filter_obj = {
            "property": prop,
            "relation": {"contains": program_id},
        }
        ids = _query_datasets(filter_obj)
        for pid in ids:
            if pid not in seen:
                collected.append(pid)
                seen.add(pid)

    if not collected:
        logger.warning(
            "[CACHE-WARM] No datasets found for program %s (properties tried: %s)",
            program_id,
            ", ".join(property_candidates),
        )

    return collected


def warm_cache(dataset_ids: Iterable[str]) -> None:
    """Warm cache for given dataset IDs."""
    cache = get_feature_cache()
    ids = list(dataset_ids)

    if not ids:
        logger.warning("[CACHE-WARM] No dataset IDs provided to warm")
        return

    start = time.time()
    successes = 0
    total_features = 0

    logger.info("[CACHE-WARM] Warming cache for %d datasets", len(ids))

    for dataset_id in tqdm(ids, desc="Warming cache", unit="dataset"):
        try:
            features = extract_dataset_features_by_type(
                dataset_page_id=dataset_id,
                use_cache=True,
                force_refresh=True,
            )
            if features:
                successes += 1
                total_features += sum(len(v) for v in features.values())

            stats = cache.get_stats()
            tqdm.write(
                f"[CACHE-WARM] Dataset {dataset_id[:8]}... | "
                f"size={stats.get('cached_datasets')} "
                f"hits={stats.get('hits')} "
                f"misses={stats.get('misses')} "
                f"evictions={stats.get('evictions')}"
            )
        except Exception as exc:
            logger.warning(
                "[CACHE-WARM] Error warming dataset %s: %r",
                dataset_id,
                exc,
            )
            continue

    elapsed = time.time() - start
    final_stats = cache.get_stats()

    print("\n=== Cache Warm Summary ===")
    print(f"Datasets requested : {len(ids)}")
    print(f"Datasets succeeded : {successes}")
    print(f"Features cached    : {total_features}")
    print(f"Time elapsed (s)   : {elapsed:.2f}")
    print(f"Cache stats        : {final_stats}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Warm the dataset feature cache for faster signature scoring",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--program-id",
        type=str,
        help="Program page ID to warm all linked datasets",
    )
    group.add_argument(
        "--all-datasets",
        action="store_true",
        help="Warm cache for all datasets in Experimental Data Assets DB",
    )
    group.add_argument(
        "--dataset-ids",
        nargs="+",
        help="Explicit dataset page IDs to warm",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.dataset_ids:
        dataset_ids = args.dataset_ids
    elif args.program_id:
        dataset_ids = fetch_dataset_ids_for_program(args.program_id)
    elif args.all_datasets:
        dataset_ids = fetch_all_dataset_ids()
    else:
        dataset_ids = []

    warm_cache(dataset_ids)


if __name__ == "__main__":
    main()

