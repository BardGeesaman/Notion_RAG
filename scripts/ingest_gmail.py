#!/usr/bin/env python3
"""
Gmail ingestion script.

Fetches emails directly from Gmail and ingests them using Postgres-only ingestion.
Replaces Zapier workflow by connecting directly to Gmail API.

Usage:
    python scripts/ingest_gmail.py                    # Ingest recent emails
    python scripts/ingest_gmail.py --all              # Ingest all inbox emails
    python scripts/ingest_gmail.py --query "is:unread"  # Custom query
    python scripts/ingest_gmail.py --days 7           # Last 7 days
"""

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from amprenta_rag.clients.gmail_client import GmailClient
from amprenta_rag.ingestion.postgres_content_ingestion import ingest_email_content
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch emails from Gmail and ingest them (Postgres-only, no Notion)"
    )
    parser.add_argument(
        "--query",
        default="in:inbox",
        help='Gmail search query (default: "in:inbox"). Examples: "is:unread", "from:sender@gmail.com"',
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=100,
        help="Maximum number of emails to fetch (default: 100)",
    )
    parser.add_argument(
        "--days",
        type=int,
        default=None,
        help="Only fetch emails from last N days (default: None = all)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Ingest all emails matching query (ignores max-results limit)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be ingested without actually ingesting",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-ingest all emails even if already ingested",
    )

    args = parser.parse_args()

    # Calculate since date if days specified
    since = None
    if args.days:
        since = datetime.now() - timedelta(days=args.days)
        logger.info("[GMAIL] Fetching emails from last %d days", args.days)

    try:
        # Initialize Gmail client
        logger.info("[GMAIL] Initializing Gmail client...")
        client = GmailClient()

        # Adjust max_results for --all flag
        max_results = None if args.all else args.max_results
        if max_results is None:
            max_results = 1000  # Large limit for --all

        # Fetch emails
        logger.info("[GMAIL] Fetching emails with query: %s", args.query)
        emails = client.fetch_emails(
            query=args.query,
            max_results=max_results,
            since=since,
        )

        if not emails:
            logger.info("[GMAIL] No emails found matching query")
            return

        logger.info("[GMAIL] Found %d email(s) to process", len(emails))

        # Process each email
        ingested_count = 0
        skipped_count = 0
        error_count = 0

        for idx, email in enumerate(emails, 1):
            email_id = email.get("id", "unknown")
            subject = email.get("subject", "(no subject)")
            from_sender = email.get("from", "")

            logger.info(
                "[GMAIL] [%d/%d] Processing: %s",
                idx,
                len(emails),
                subject[:60],
            )

            if args.dry_run:
                logger.info(
                    "[GMAIL] [DRY RUN] Would ingest: %s (from: %s)",
                    subject,
                    from_sender,
                )
                continue

            try:
                # Extract email body
                body = email.get("body", "")
                if not body or len(body.strip()) < 50:
                    logger.debug(
                        "[GMAIL] Skipping email %s: body too short",
                        email_id[:8],
                    )
                    skipped_count += 1
                    continue

                # Extract tags from labels (Gmail labels)
                labels = email.get("labels", [])
                tags = [
                    label
                    for label in labels
                    if label not in ["INBOX", "UNREAD", "IMPORTANT", "SENT"]
                ]

                # Extract additional metadata
                metadata = {
                    "to": email.get("to", ""),
                    "cc": email.get("cc", ""),
                    "date": email.get("date", ""),
                    "gmail_labels": labels,
                    "thread_id": email.get("thread_id", ""),
                }

                # Ingest email using Postgres-only ingestion (idempotent)
                embedding_ids = ingest_email_content(
                    email_content=body,
                    title=subject,
                    from_sender=from_sender,
                    email_id=f"gmail_{email_id}",
                    tags=tags,
                    metadata=metadata,
                    force=args.force,
                )

                if embedding_ids:
                    ingested_count += 1
                    logger.info(
                        "[GMAIL] ✅ Ingested email: %s (%d chunks)",
                        subject[:60],
                        len(embedding_ids),
                    )
                else:
                    skipped_count += 1
                    logger.info(
                        "[GMAIL] ⏭️  Skipped email: %s (already ingested)",
                        subject[:60],
                    )

            except Exception as e:
                error_count += 1
                logger.error(
                    "[GMAIL] ❌ Error ingesting email %s: %r",
                    email_id[:8],
                    e,
                )

        # Summary
        logger.info("[GMAIL] " + "=" * 60)
        logger.info("[GMAIL] Ingestion complete!")
        logger.info("[GMAIL]   Total emails: %d", len(emails))
        if not args.dry_run:
            logger.info("[GMAIL]   Ingested: %d", ingested_count)
            logger.info("[GMAIL]   Skipped: %d", skipped_count)
            logger.info("[GMAIL]   Errors: %d", error_count)

    except Exception as e:
        logger.error("[GMAIL] Fatal error: %r", e)
        sys.exit(1)


if __name__ == "__main__":
    main()

