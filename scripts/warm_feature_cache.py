#!/usr/bin/env python3
"""
Warm the dataset feature cache for faster signature scoring.

Supports warming by:
- Program ID (Postgres)
- All datasets (Postgres)
- Explicit dataset IDs

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Set
from uuid import UUID

from tqdm import tqdm

# Ensure repository imports work when run as a script
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Program
from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _query_datasets() -> List[str]:
    """Query Postgres datasets and return dataset IDs."""
    with db_session() as db:
        query = db.query(Dataset)
        datasets = query.limit(500).all()
        return [str(d.id) for d in datasets]


def fetch_all_dataset_ids() -> List[str]:
    """List all dataset IDs from Postgres."""
    datasets = _query_datasets()
    if not datasets:
        logger.warning("[CACHE-WARM] No datasets found in Postgres")
    return datasets


def fetch_dataset_ids_for_program(program_id: str) -> List[str]:
    """List dataset IDs linked to a Program in Postgres."""
    with db_session() as db:
        # Try to parse as UUID
        try:
            prog_uuid = UUID(program_id)
            program = db.query(Program).filter(Program.id == prog_uuid).first()
        except ValueError:
            # Not a valid UUID, try by name
            program = db.query(Program).filter(Program.name.ilike(f"%{program_id}%")).first()
        
        if not program:
            logger.warning("[CACHE-WARM] Program not found: %s", program_id)
            return []
        
        # Get datasets linked to this program
        datasets = program.datasets if hasattr(program, 'datasets') else []
        return [str(d.id) for d in datasets]


def warm_cache_for_datasets(dataset_ids: List[str], batch_size: int = 50) -> None:
    """
    Warm the feature cache for a list of datasets.

    Args:
        dataset_ids: List of dataset IDs (UUIDs as strings)
        batch_size: Number of datasets to process at once
    """
    cache = get_feature_cache()
    total = len(dataset_ids)

    logger.info("[CACHE-WARM] Warming cache for %d datasets...", total)

    warmed = 0
    skipped = 0
    errors = 0

    for dataset_id in tqdm(dataset_ids, desc="Warming cache"):
        # Check if already cached
        if cache.get(dataset_id):
            skipped += 1
            continue

        try:
            features = extract_dataset_features_by_type(dataset_id)
            if features:
                cache.put(dataset_id, features)
                warmed += 1
            else:
                skipped += 1
        except Exception as exc:
            logger.debug("[CACHE-WARM] Error extracting features for %s: %r", dataset_id, exc)
            errors += 1

    logger.info(
        "[CACHE-WARM] Complete. Warmed: %d, Skipped: %d, Errors: %d",
        warmed,
        skipped,
        errors,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Warm the dataset feature cache for faster signature scoring."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--program",
        metavar="PROGRAM_ID",
        help="Warm cache for datasets linked to a program",
    )
    group.add_argument(
        "--all-datasets",
        action="store_true",
        help="Warm cache for all datasets",
    )
    group.add_argument(
        "--dataset-ids",
        nargs="+",
        metavar="ID",
        help="Warm cache for specific dataset IDs",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Number of datasets to process at once (default: 50)",
    )
    args = parser.parse_args()

    if args.program:
        dataset_ids = fetch_dataset_ids_for_program(args.program)
        if not dataset_ids:
            print(f"No datasets found for program: {args.program}")
            sys.exit(1)
    elif args.all_datasets:
        dataset_ids = fetch_all_dataset_ids()
        if not dataset_ids:
            print("No datasets found in database")
            sys.exit(1)
    else:
        dataset_ids = args.dataset_ids

    print(f"Found {len(dataset_ids)} datasets to warm")
    warm_cache_for_datasets(dataset_ids, batch_size=args.batch_size)


if __name__ == "__main__":
    main()
