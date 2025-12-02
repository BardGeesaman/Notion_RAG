#!/usr/bin/env python3
# scripts/update_collection_universe.py

import argparse
from amprenta_rag.maintenance import update_collection_universe


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Incrementally update and delete-sync a Zotero collection universe."
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key (e.g. 3RGXZTAY).",
    )
    args = parser.parse_args()
    update_collection_universe(args.collection_key)


if __name__ == "__main__":
    main()