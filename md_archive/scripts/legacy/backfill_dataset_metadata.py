#!/usr/bin/env python3
"""
Backfill script to extract and update scientific metadata for existing Experimental Data Assets pages.

This script:
1. Finds all pages in the Experimental Data Assets database
2. Filters for pages with mwTab data but missing metadata (empty Summary or Disease)
3. Extracts metadata from mwTab content
4. Updates the Notion pages with extracted metadata

Usage:
    python scripts/backfill_dataset_metadata.py
    python scripts/backfill_dataset_metadata.py --dry-run
    python scripts/backfill_dataset_metadata.py --page-id <page_id>
"""

import argparse
import sys
from pathlib import Path

import requests

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.dataset_notion_utils import (
    update_dataset_scientific_metadata)
from amprenta_rag.ingestion.mwtab_extraction import (
    extract_metadata_from_mwtab, extract_mwtab_from_page_content)
from amprenta_rag.ingestion.notion_pages import extract_page_content
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def fetch_all_dataset_pages() -> list:
    """Fetch all pages from Experimental Data Assets database."""
    cfg = get_config().notion
    if not cfg.exp_data_db_id:
        logger.error("Experimental Data Assets DB ID not configured")
        return []

    url = f"{cfg.base_url}/databases/{cfg.exp_data_db_id}/query"
    all_pages = []
    next_cursor = None

    while True:
        payload = {"page_size": 100}
        if next_cursor:
            payload["start_cursor"] = next_cursor

        try:
            resp = requests.post(
                url, headers=notion_headers(), json=payload, timeout=30
            )
            resp.raise_for_status()
            data = resp.json()
            pages = data.get("results", [])
            all_pages.extend(pages)

            if not data.get("has_more"):
                break
            next_cursor = data.get("next_cursor")
        except Exception as e:
            logger.error("Error fetching dataset pages: %r", e)
            break

    return all_pages


def page_needs_metadata(page: dict) -> bool:
    """
    Check if page needs metadata backfill.

    Returns True if page is missing Summary or Disease fields.
    Note: We'll check for mwTab data later during processing.
    """
    props = page.get("properties", {}) or {}

    # Check if Summary is empty or only has MW Study ID
    summary_prop = props.get("Summary", {})
    summary_rich = summary_prop.get("rich_text", []) or []
    summary_text = "".join(rt.get("plain_text", "") for rt in summary_rich).strip()
    # Summary is considered empty if it only contains "MW Study ID: ..." or is empty
    has_summary = bool(summary_text and not summary_text.startswith("MW Study ID:"))

    # Check if Disease is empty
    disease_prop = props.get("Disease", {})
    disease_multi = disease_prop.get("multi_select", []) or []
    has_disease = len(disease_multi) > 0

    # Page needs metadata if it's missing Summary or Disease
    return not (has_summary and has_disease)


def process_page(page_id: str, dry_run: bool = False) -> bool:
    """
    Process a single dataset page: extract metadata and update if needed.

    Returns:
        True if metadata was updated (or would be updated in dry-run), False otherwise
    """
    logger.info("Processing page %s", page_id)

    try:
        # Extract page content
        page_id_clean = page_id.replace("-", "")
        page_content = extract_page_content(page_id_clean)

        if not page_content:
            logger.warning("No content found for page %s", page_id)
            return False

        # Extract mwTab data
        mwtab_data = extract_mwtab_from_page_content(page_content)
        if not mwtab_data:
            logger.info("No mwTab data found for page %s", page_id)
            return False

        # Extract metadata
        metadata = extract_metadata_from_mwtab(mwtab_data)

        # Check if we have any metadata to set
        has_metadata = any(
            [
                metadata.get("model_systems"),
                metadata.get("disease_terms"),
                metadata.get("matrix_terms"),
                metadata.get("methods"),
                metadata.get("summary"),
                metadata.get("results"),
            ]
        )

        if not has_metadata:
            logger.info("No extractable metadata for page %s", page_id)
            return False

        if dry_run:
            logger.info("[DRY-RUN] Would update page %s with metadata:", page_id)
            logger.info("  Model Systems: %r", metadata.get("model_systems"))
            logger.info("  Disease: %r", metadata.get("disease_terms"))
            logger.info("  Matrix: %r", metadata.get("matrix_terms"))
            logger.info("  Methods: %s", metadata.get("methods", "")[:100])
            logger.info("  Summary: %s", metadata.get("summary", "")[:100])
            return True
        else:
            # Update page
            update_dataset_scientific_metadata(
                page_id=page_id,
                metadata=metadata,
            )
            return True

    except Exception as e:
        logger.error("Error processing page %s: %r", page_id, e)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Backfill scientific metadata for Experimental Data Assets pages"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be updated without making changes",
    )
    parser.add_argument(
        "--page-id",
        help="Process only this specific page ID (with or without dashes)",
    )

    args = parser.parse_args()

    if args.page_id:
        # Process single page
        page_id = args.page_id
        if "-" not in page_id:
            # Add dashes if missing
            page_id = f"{page_id[:8]}-{page_id[8:12]}-{page_id[12:16]}-{page_id[16:20]}-{page_id[20:]}"

        logger.info("Processing single page: %s", page_id)
        process_page(page_id, dry_run=args.dry_run)
    else:
        # Fetch all pages and filter
        logger.info("Fetching all Experimental Data Assets pages...")
        all_pages = fetch_all_dataset_pages()
        logger.info("Found %d pages total", len(all_pages))

        # Filter pages that need metadata
        pages_to_process = [p for p in all_pages if page_needs_metadata(p)]
        logger.info("Found %d pages that need metadata backfill", len(pages_to_process))

        if args.dry_run:
            logger.info("[DRY-RUN] Would process %d pages", len(pages_to_process))

        # Process each page
        updated_count = 0
        for page in pages_to_process:
            page_id = page.get("id", "")
            if process_page(page_id, dry_run=args.dry_run):
                updated_count += 1

        logger.info("Processed %d pages", updated_count)


if __name__ == "__main__":
    main()
