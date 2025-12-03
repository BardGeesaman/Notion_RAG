#!/usr/bin/env python3
# scripts/ingest_email.py

import argparse

from amprenta_rag.ingestion import batch_ingest_emails, delete_email_and_chunks


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-ingest emails from Email & Notes Inbox into RAG Engine + Pinecone."
    )
    parser.add_argument(
        "--parent-type",
        default="Email",
        help='Parent type label for Notion (default: "Email").',
    )
    parser.add_argument(
        "--delete-email",
        help=(
            "Delete an email and its chunks by page ID "
            "(removes from Notion and Pinecone)."
        ),
    )

    args = parser.parse_args()

    if args.delete_email:
        delete_email_and_chunks(args.delete_email)
    else:
        batch_ingest_emails(parent_type=args.parent_type)


if __name__ == "__main__":
    main()
