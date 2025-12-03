#!/usr/bin/env python3
# scripts/ingest_collection.py

import argparse

from amprenta_rag.ingestion import incremental_ingest_collection


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Incrementally ingest a Zotero collection into Notion + Pinecone."
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

    incremental_ingest_collection(
        collection_key=args.collection_key,
        parent_type=args.parent_type,
    )


if __name__ == "__main__":
    main()
