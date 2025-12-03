#!/usr/bin/env python3
"""
CLI for ingesting internal lipidomics datasets.

Usage:
    python scripts/ingest_lipidomics.py --file <path> --create-page
    python scripts/ingest_lipidomics.py --file <path> --dataset-page-id <id>
    python scripts/ingest_lipidomics.py --file <path> --create-page --program-id <id> --experiment-id <id>
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest an internal lipidomics dataset file into Notion and RAG."
    )
    parser.add_argument(
        "--file",
        required=True,
        help="Path to the lipidomics file (CSV/TSV)",
    )
    parser.add_argument(
        "--dataset-page-id",
        default=None,
        help="Existing Experimental Data Assets page ID to attach to",
    )
    parser.add_argument(
        "--create-page",
        action="store_true",
        help="Create a new Experimental Data Asset page (required if --dataset-page-id not provided)",
    )
    parser.add_argument(
        "--program-id",
        dest="program_ids",
        action="append",
        default=[],
        help="Program page ID to link (can be specified multiple times)",
    )
    parser.add_argument(
        "--experiment-id",
        dest="experiment_ids",
        action="append",
        default=[],
        help="Experiment page ID to link (can be specified multiple times)",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.dataset_page_id and not args.create_page:
        parser.error(
            "Either --dataset-page-id must be provided or --create-page must be specified"
        )

    try:
        # Ingest the file
        page_id = ingest_lipidomics_file(
            file_path=args.file,
            notion_page_id=args.dataset_page_id,
            create_page=args.create_page,
            program_ids=args.program_ids if args.program_ids else None,
            experiment_ids=args.experiment_ids if args.experiment_ids else None,
        )

        print(f"\n✅ Ingestion complete!")
        print(f"   Dataset page ID: {page_id}")
        print(f"   File: {args.file}")
        if args.program_ids:
            print(f"   Linked to {len(args.program_ids)} program(s)")
        if args.experiment_ids:
            print(f"   Linked to {len(args.experiment_ids)} experiment(s)")

    except Exception as e:
        logger.error("[INGEST][LIPIDOMICS] Ingestion failed: %r", e)
        print(f"\n❌ Ingestion failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

