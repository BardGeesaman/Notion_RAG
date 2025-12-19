#!/usr/bin/env python3
"""
CLI for ingesting Experiment ELN pages into RAG.

Usage:
    python scripts/ingest_experiment.py --experiment-page-id <id>
    python scripts/ingest_experiment.py --experiment-page-id <id> --force
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.experiments_ingestion import ingest_experiment
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Ingest a single Experiment ELN page from Notion into RAG."
    )
    parser.add_argument(
        "--experiment-page-id",
        required=True,
        help="Notion page ID (with or without dashes) for the Experiment page.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-ingestion (bypass caching, for future use)",
    )
    args = parser.parse_args()

    try:
        ingest_experiment(
            exp_page_id=args.experiment_page_id,
            parent_type="Experiment",
            force=args.force,
        )

        print("\n✅ Ingestion complete!")
        print(f"   Experiment page ID: {args.experiment_page_id}")
        print("   Check logs above for chunk count and embedding details")

    except Exception as e:
        logger.error("[INGEST][EXPERIMENT] Ingestion failed: %r", e)
        print(f"\n❌ Ingestion failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
