#!/usr/bin/env python3
"""
Zotero ingestion script (Postgres-only).

Fetches items from Zotero collections and ingests them using Postgres-only ingestion.
No Notion dependencies - direct to Pinecone.

Usage:
    python scripts/ingest_zotero_postgres.py --collection-key COLLECTION_KEY
    python scripts/ingest_zotero_postgres.py --collection-key COLLECTION_KEY --max-items 5
    python scripts/ingest_zotero_postgres.py --collection-key COLLECTION_KEY --dry-run
"""

import argparse
import sys
from pathlib import Path

from tqdm import tqdm
# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from amprenta_rag.ingestion.postgres_content_ingestion import ingest_literature_content
from amprenta_rag.ingestion.text_extraction import extract_text_from_attachment_bytes
from amprenta_rag.ingestion.zotero_collection import fetch_collection_items
from amprenta_rag.ingestion.zotero_api import (
    ZoteroItem,
    download_zotero_file,
    fetch_zotero_attachments,
    fetch_zotero_item,
    fetch_zotero_notes,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_authors_from_item(data: dict) -> list[str]:
    """Extract authors from Zotero item data."""
    creators = data.get("creators", [])
    authors = []
    for creator in creators:
        if creator.get("creatorType") == "author":
            first = creator.get("firstName", "")
            last = creator.get("lastName", "")
            if first or last:
                authors.append(f"{first} {last}".strip())
    return authors


def process_zotero_item(
    item_key: str, dry_run: bool = False, force: bool = False
) -> tuple[int, int, int]:
    """
    Process a single Zotero item: fetch metadata, attachments, notes, and ingest.

    Returns:
        Tuple of (ingested_count, skipped_count, error_count)
    """
    ingested_count = 0
    skipped_count = 0
    error_count = 0

    try:
        # Fetch item metadata
        item: ZoteroItem = fetch_zotero_item(item_key)
        logger.info(
            "[ZOTERO] Processing item: %s (%s)",
            item.title[:60],
            item.item_type,
        )

        if dry_run:
            logger.info(
                "[ZOTERO] [DRY RUN] Would ingest: %s (key: %s)",
                item.title,
                item_key,
            )
            return 0, 0, 0

        # Fetch full item data for authors
        from amprenta_rag.ingestion.zotero_api import _zotero_base, _zotero_headers
        import requests

        url = f"{_zotero_base()}/items/{item_key}?include=data"
        resp = requests.get(url, headers=_zotero_headers(), timeout=30)
        resp.raise_for_status()
        item_data = resp.json().get("data", {})

        authors = extract_authors_from_item(item_data)

        # Build content from abstract, notes, and attachments
        content_parts = []

        # Add abstract if available
        if item.abstract and len(item.abstract.strip()) > 50:
            content_parts.append(f"Abstract:\n{item.abstract}\n")

        # Process notes
        try:
            notes = fetch_zotero_notes(item_key)
            for note in notes:
                note_content = note.get("note", "")
                if note_content:
                    # Remove HTML tags (simple strip)
                    import re

                    text_content = re.sub(r"<[^>]+>", "", note_content)
                    if text_content.strip():
                        content_parts.append(f"Note:\n{text_content}\n")
        except Exception as e:
            logger.debug("[ZOTERO] Error fetching notes: %r", e)

        # Process attachments (PDFs, etc.)
        try:
            attachments = fetch_zotero_attachments(item_key)
            for att in attachments:
                att_key = att.get("key")
                filename = att.get("filename", "")
                content_type = att.get("contentType", "")

                # Skip HTML snapshots
                if "snapshot" in filename.lower() or "html" in content_type.lower():
                    logger.debug("[ZOTERO] Skipping snapshot/HTML attachment: %s", filename)
                    continue

                # Skip non-readable attachments
                if not filename or (
                    not filename.endswith(".pdf")
                    and not content_type.startswith("text/")
                    and "pdf" not in content_type.lower()
                ):
                    logger.debug(
                        "[ZOTERO] Skipping non-text attachment: %s (%s)",
                        filename,
                        content_type,
                    )
                    continue

                try:
                    logger.info(
                        "[ZOTERO] Processing attachment: %s",
                        filename or att_key[:8],
                    )
                    file_bytes = download_zotero_file(att_key)
                    text = extract_text_from_attachment_bytes(
                        content_type, file_bytes, filename=filename
                    )

                    if text and len(text.strip()) > 100:
                        content_parts.append(f"Attachment ({filename}):\n{text}\n")
                except Exception as e:
                    logger.warning(
                        "[ZOTERO] Error processing attachment %s: %r",
                        filename or att_key[:8],
                        e,
                    )
                    skipped_count += 1
        except Exception as e:
            logger.debug("[ZOTERO] Error fetching attachments: %r", e)

        # Combine all content
        full_content = "\n".join(content_parts)

        if not full_content or len(full_content.strip()) < 100:
            logger.warning(
                "[ZOTERO] Item %s has insufficient content (%d chars), skipping",
                item_key,
                len(full_content) if full_content else 0,
            )
            skipped_count += 1
            return ingested_count, skipped_count, error_count

        # Build metadata
        metadata = {
            "item_type": item.item_type,
            "journal": item.journal or "",
            "date": item.date or "",
            "url": item.url or "",
            "tags": item.tags,
        }

        # Ingest using Postgres-only ingestion (idempotent)
        embedding_ids = ingest_literature_content(
            literature_content=full_content,
            title=item.title,
            authors=authors,
            doi=item.doi,
            zotero_key=item_key,
            metadata=metadata,
            force=force,
        )

        if embedding_ids:
            ingested_count = 1
            logger.info(
                "[ZOTERO] ✅ Ingested item: %s (%d chunks)",
                item.title[:60],
                len(embedding_ids),
            )
        else:
            skipped_count = 1
            logger.info(
                "[ZOTERO] ⏭️  Skipped item: %s (already ingested)",
                item.title[:60],
            )

    except Exception as e:
        error_count += 1
        logger.error("[ZOTERO] ❌ Error processing item %s: %r", item_key[:8], e)

    return ingested_count, skipped_count, error_count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch items from Zotero collection and ingest (Postgres-only, no Notion)"
    )
    parser.add_argument(
        "--collection-key",
        required=True,
        help="Zotero collection key (e.g. 3RGXZTAY)",
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=None,
        help="Maximum number of items to process (default: all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be ingested without actually ingesting",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-ingest all items even if already ingested",
    )

    args = parser.parse_args()

    try:
        # Fetch collection items
        logger.info(
            "[ZOTERO] Fetching items from collection: %s",
            args.collection_key,
        )
        items = fetch_collection_items(args.collection_key)

        if not items:
            logger.info("[ZOTERO] No items found in collection")
            return

        # Filter to literature-like items
        allowed_types = {
            "journalArticle",
            "book",
            "bookSection",
            "report",
            "preprint",
            "conferencePaper",
        }
        literature_items = [
            item
            for item in items
            if item.get("data", {}).get("itemType") in allowed_types
        ]

        if not literature_items:
            logger.info("[ZOTERO] No literature items found in collection")
            return

        # Limit items if requested
        if args.max_items:
            literature_items = literature_items[: args.max_items]

        logger.info(
            "[ZOTERO] Found %d literature item(s) to process",
            len(literature_items),
        )

        # Process each item
        total_ingested = 0
        total_skipped = 0
        total_errors = 0

        for idx, item in enumerate(tqdm(literature_items, desc="Ingesting Zotero items"), 1):
            item_key = item.get("data", {}).get("key")
            if not item_key:
                continue

            item_title = item.get("data", {}).get("title", "(untitled)")[:60]
            logger.info(
                "[ZOTERO] [%d/%d] Processing: %s",
                idx,
                len(literature_items),
                item_title,
            )

            ingested, skipped, errors = process_zotero_item(
                item_key, dry_run=args.dry_run, force=args.force
            )
            total_ingested += ingested
            total_skipped += skipped
            total_errors += errors

        # Summary
        logger.info("[ZOTERO] " + "=" * 60)
        logger.info("[ZOTERO] Ingestion complete!")
        logger.info("[ZOTERO]   Total items: %d", len(literature_items))
        if not args.dry_run:
            logger.info("[ZOTERO]   Ingested: %d", total_ingested)
            logger.info("[ZOTERO]   Skipped: %d", total_skipped)
            logger.info("[ZOTERO]   Errors: %d", total_errors)

    except Exception as e:
        logger.error("[ZOTERO] Fatal error: %r", e)
        sys.exit(1)


if __name__ == "__main__":
    main()

