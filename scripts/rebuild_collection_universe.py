#!/usr/bin/env python3
# scripts/rebuild_collection_universe.py

import argparse

from amprenta_rag.maintenance import rebuild_collection_universe


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rebuild the entire RAG universe for a Zotero collection (full resync + cleanup)."
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key (e.g. 3RGXZTAY).",
    )
    args = parser.parse_args()
    rebuild_collection_universe(args.collection_key)


if __name__ == "__main__":
    main()
