#!/usr/bin/env python3
# scripts/resync_collection.py

import argparse

from amprenta_rag.ingestion import resync_collection


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fully resync a Zotero collection (delete + re-ingest)."
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key (e.g. 3RGXZTAY).",
    )
    parser.add_argument(
        "--parent-type",
        default="Literature",
        help='Parent type label for Notion (default: "Literature").',
    )
    args = parser.parse_args()

    resync_collection(
        collection_key=args.collection_key,
        parent_type=args.parent_type,
    )


if __name__ == "__main__":
    main()
