#!/usr/bin/env python3
"""
Manage the dataset feature cache: stats, clear, invalidate, export, import.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, Iterable

# Ensure repository imports resolve when invoked as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _format_stats(stats: Dict[str, int]) -> str:
    hits = stats.get("hits", 0)
    misses = stats.get("misses", 0)
    total = hits + misses
    hit_rate = (hits / total) if total else 0.0

    rows = [
        ("cached_datasets", stats.get("cached_datasets", 0)),
        ("total_features", stats.get("total_features", 0)),
        ("hits", hits),
        ("misses", misses),
        ("evictions", stats.get("evictions", 0)),
        ("hit_rate", f"{hit_rate:.1%}"),
    ]

    width = max(len(k) for k, _ in rows) + 2
    lines = ["Cache Statistics:"]
    for key, val in rows:
        lines.append(f"  {key.ljust(width)}{val}")
    return "\n".join(lines)


def cmd_stats() -> None:
    cache = get_feature_cache()
    stats = cache.get_stats()
    print(_format_stats(stats))


def cmd_clear(force: bool) -> None:
    if not force:
        print("Refusing to clear cache without --force")
        return
    cache = get_feature_cache()
    cache.clear()
    logger.info("[CACHE-MANAGE] Cleared entire cache")


def cmd_invalidate(dataset_id: str) -> None:
    cache = get_feature_cache()
    cache.invalidate(dataset_id)
    logger.info("[CACHE-MANAGE] Invalidated dataset %s", dataset_id)


def cmd_export(file_path: str) -> None:
    cache = get_feature_cache()

    # Snapshot under lock to avoid race
    with cache._lock:  # type: ignore[attr-defined]
        entries = {}
        for ds_id, entry in cache.cache.items():
            entries[ds_id] = {
                "features": {k: list(v) for k, v in entry["features"].items()},
                "timestamp": entry.get("timestamp", time.time()),
                "ttl": entry.get("ttl", cache.default_ttl),
            }

    payload = {
        "exported_at": time.time(),
        "entries": entries,
    }

    try:
        with open(file_path, "w") as f:
            json.dump(payload, f, indent=2)
        logger.info(
            "[CACHE-MANAGE] Exported %d entries to %s",
            len(entries),
            file_path,
        )
    except Exception as exc:
        logger.warning("[CACHE-MANAGE] Failed to export cache: %r", exc)


def cmd_import(file_path: str) -> None:
    cache = get_feature_cache()
    try:
        with open(file_path, "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        logger.warning("[CACHE-MANAGE] Import file not found: %s", file_path)
        return
    except Exception as exc:
        logger.warning("[CACHE-MANAGE] Failed to read import file: %r", exc)
        return

    entries = data.get("entries", {}) or {}
    loaded = 0

    with cache._lock:  # type: ignore[attr-defined]
        for ds_id, entry in entries.items():
            features_raw = entry.get("features", {}) or {}
            features = {k: set(v or []) for k, v in features_raw.items()}
            timestamp = entry.get("timestamp", time.time())
            ttl = entry.get("ttl", cache.default_ttl)

            cache.cache[ds_id] = {
                "features": features,
                "timestamp": timestamp,
                "ttl": ttl,
            }
            cache.cache.move_to_end(ds_id)
            cache._evict_if_needed_unlocked()  # type: ignore[attr-defined]
            loaded += 1

    logger.info(
        "[CACHE-MANAGE] Imported %d entries from %s",
        loaded,
        file_path,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Manage the dataset feature cache",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--stats", action="store_true", help="Show cache statistics")
    group.add_argument(
        "--clear",
        action="store_true",
        help="Clear the entire cache (requires --force)",
    )
    group.add_argument(
        "--invalidate",
        metavar="DATASET_ID",
        help="Invalidate a specific dataset ID",
    )
    group.add_argument(
        "--export",
        metavar="FILE",
        help="Export cache to JSON file",
    )
    group.add_argument(
        "--import",
        dest="import_file",
        metavar="FILE",
        help="Import cache from JSON file",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Confirm destructive operations (used with --clear)",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.stats:
        cmd_stats()
    elif args.clear:
        cmd_clear(force=args.force)
    elif args.invalidate:
        cmd_invalidate(args.invalidate)
    elif args.export:
        cmd_export(args.export)
    elif args.import_file:
        cmd_import(args.import_file)
    else:
        raise SystemExit("No command provided")


if __name__ == "__main__":
    main()

