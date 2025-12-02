#!/usr/bin/env python3
# scripts/rag_reindex_all.py

import argparse

from amprenta_rag.maintenance import rebuild_collection_universe
from amprenta_rag.ingestion import batch_ingest_emails


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reindex the entire RAG universe (Zotero collection + Emails)."
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key to rebuild (e.g. 3RGXZTAY).",
    )
    parser.add_argument(
        "--skip-emails",
        action="store_true",
        help="If set, only rebuild Zotero side; skip email ingestion.",
    )
    args = parser.parse_args()

    print("\n🚀 Rebuilding Zotero side...")
    rebuild_collection_universe(args.collection_key)

    if not args.skip_emails:
        print("\n📨 Ingesting emails/notes from Email & Notes Inbox...")
        batch_ingest_emails(parent_type="Email")

    print("\n🎉 Full RAG reindex complete.\n")


if __name__ == "__main__":
    main()