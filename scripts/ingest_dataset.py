#!/usr/bin/env python3

import argparse

from amprenta_rag.ingestion.dataset_ingestion import ingest_dataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-page-id", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    ingest_dataset(args.dataset_page_id, force=args.force)


if __name__ == "__main__":
    main()
