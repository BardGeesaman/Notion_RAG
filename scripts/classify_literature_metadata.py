#!/usr/bin/env python3
# scripts/classify_literature_metadata.py

import argparse

from amprenta_rag.metadata import classify_and_update_all_literature


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify semantic + lipid metadata for Literature DB pages and update Notion."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="If set, only classify up to this many pages (for testing).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="If set, reclassify and update pages even if they already have semantic metadata.",
    )
    args = parser.parse_args()

    classify_and_update_all_literature(limit=args.limit, force=args.force)


if __name__ == "__main__":
    main()