#!/usr/bin/env python

import argparse

from amprenta_rag.ingestion.experiments_ingestion import ingest_experiment


def main():
    parser = argparse.ArgumentParser(
        description="Ingest a single Experiment page from Notion into Pinecone."
    )
    parser.add_argument(
        "--experiment-page-id",
        required=True,
        help="Notion page ID (with or without dashes) for the Experiment page.",
    )
    args = parser.parse_args()

    ingest_experiment(args.experiment_page_id, parent_type="Experiment")


if __name__ == "__main__":
    main()
