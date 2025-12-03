#!/usr/bin/env python3
# scripts/sync_collection_state.py

import argparse

from amprenta_rag.maintenance import sync_collection_state


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sync collection-defined universe: keep only items in given Zotero collection."
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key (e.g. 3RGXZTAY).",
    )
    args = parser.parse_args()
    sync_collection_state(args.collection_key)


if __name__ == "__main__":
    main()
