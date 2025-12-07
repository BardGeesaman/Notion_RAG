"""
Main Zotero ingestion orchestration.

This module provides the main ingestion function that orchestrates the complete
Zotero item ingestion pipeline: fetching items, processing attachments/notes,
feature extraction, signature detection, and Notion updates.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import List

from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text,
    link_features_to_notion_items,
)
from amprenta_rag.ingestion.metadata_semantic import get_literature_semantic_metadata
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.ingestion.notion_pages import (
#     ensure_literature_page,
#     update_literature_page,
# )

def ensure_literature_page(item, parent_type: str) -> str:
    """DEPRECATED: Notion support removed. Returns placeholder page ID."""
    logger.debug("[ZOTERO][INGESTION] ensure_literature_page() deprecated - Notion support removed")
    return "placeholder-page-id"

def update_literature_page(page_id: str, when_iso: str) -> None:
    """DEPRECATED: Notion support removed. Does nothing."""
    logger.debug("[ZOTERO][INGESTION] update_literature_page() deprecated - Notion support removed")
    return
from amprenta_rag.ingestion.signature_integration import (
    detect_and_ingest_signatures_from_content,
)
from amprenta_rag.ingestion.zotero.attachment_processing import process_attachments
from amprenta_rag.ingestion.zotero.note_processing import process_notes
from amprenta_rag.ingestion.zotero_api import (
    ZoteroItem,
    fetch_zotero_attachments,
    fetch_zotero_item,
    fetch_zotero_notes,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def ingest_zotero_item(
    item_key: str,
    parent_type: str = "Literature",
    force: bool = False,
) -> None:
    """
    Ingest a single Zotero item into Notion + Pinecone, with
    attachment-level and note-level idempotency.

    For each attachment OR note on the Zotero item:
      - Skip if already ingested with same content (md5 / hash) unless force=True
      - Attachments: download file, extract text (PDF or text/*)
      - Notes: HTML → text
      - Chunk + embed
      - Create RAG chunk pages in Notion
      - Upsert vectors with rich metadata

    If all attachments/notes are up to date, this is a no-op.
    """
    logger.info("[INGEST][ZOTERO] Ingesting Zotero item %s", item_key)

    # 1) Item metadata
    try:
        item: ZoteroItem = fetch_zotero_item(item_key)
    except Exception as e:
        logger.error(
            "[INGEST][ZOTERO] Zotero API error fetching item %s: %r", item_key, e
        )
        raise

    # 2) Ensure Literature page exists
    try:
        parent_page_id = ensure_literature_page(item, parent_type)
    except Exception as e:
        logger.error(
            "[INGEST][ZOTERO] Notion API error ensuring literature page for item %s: %r",
            item_key,
            e,
        )
        raise

    # 2b) Read semantic + lipid metadata from the Literature page
    try:
        base_meta = get_literature_semantic_metadata(parent_page_id, item)
    except Exception as e:
        logger.error(
            "[INGEST][ZOTERO] Error reading semantic metadata for page %s: %r",
            parent_page_id,
            e,
        )
        raise

    # 3) Fetch children
    try:
        attachments = fetch_zotero_attachments(item_key)
        notes = fetch_zotero_notes(item_key)
    except Exception as e:
        logger.error(
            "[INGEST][ZOTERO] Zotero API error fetching children for item %s: %r",
            item_key,
            e,
        )
        raise

    if not attachments and not notes:
        logger.info(
            "[INGEST][ZOTERO] No attachments or notes for %s; nothing to ingest.",
            item_key,
        )
        return

    logger.info(
        "[INGEST][ZOTERO] Found %d attachment(s) and %d note(s) for Zotero item %s",
        len(attachments),
        len(notes),
        item_key,
    )

    now = datetime.now(timezone.utc).isoformat()
    all_text_parts: List[str] = []

    # Process attachments
    att_texts, att_ingested = process_attachments(
        item=item,
        attachments=attachments,
        parent_page_id=parent_page_id,
        parent_type=parent_type,
        base_meta=base_meta,
        force=force,
    )
    all_text_parts.extend(att_texts)

    # Process notes
    note_texts, note_ingested = process_notes(
        item=item,
        notes=notes,
        parent_page_id=parent_page_id,
        parent_type=parent_type,
        base_meta=base_meta,
        force=force,
    )
    all_text_parts.extend(note_texts)

    any_ingested = att_ingested or note_ingested

    if any_ingested:
        logger.info(
            "[INGEST][ZOTERO] Updating Literature page status in Notion for item %s",
            item_key,
        )
        try:
            update_literature_page(parent_page_id, now)
        except Exception as e:
            logger.error(
                "[INGEST][ZOTERO] Notion API error updating literature page %s: %r",
                parent_page_id,
                e,
            )
            raise

        # Extract and link metabolite features from all ingested content
        try:
            if all_text_parts:
                combined_text = "\n\n".join(all_text_parts)
                feature_names = extract_features_from_text(combined_text)
                if feature_names:
                    logger.info(
                        "[INGEST][LITERATURE] Extracted %d metabolite feature(s) from item %s",
                        len(feature_names),
                        item_key,
                    )
                    link_features_to_notion_items(
                        feature_names=feature_names,
                        item_page_id=parent_page_id,
                        item_type="literature",
                    )
                else:
                    logger.debug(
                        "[INGEST][LITERATURE] No metabolite features found in item %s",
                        item_key,
                    )
        except Exception as e:
            # Log but don't raise - feature extraction is non-critical
            logger.warning(
                "[INGEST][LITERATURE] Error extracting/linking features for item %s: %r",
                item_key,
                e,
            )

        # Detect and ingest signatures from literature content
        try:
            if all_text_parts:
                combined_text = "\n\n".join(all_text_parts)

                # Get source metadata for signature inference
                base_meta = get_literature_semantic_metadata(parent_page_id, item)
                source_metadata = {
                    "diseases": base_meta.get("diseases", []),
                    "matrix": base_meta.get("matrix", []),
                    "model_systems": base_meta.get("model_systems", []),
                }

                # No attached files for literature (text extracted, files not stored locally)
                attachment_paths: List[Path] = []

                detect_and_ingest_signatures_from_content(
                    all_text_content=combined_text,
                    attachment_paths=attachment_paths,
                    source_page_id=parent_page_id,
                    source_type="literature",
                    source_metadata=source_metadata,
                    source_name=item.title or f"Zotero Item {item_key}",
                )
        except Exception as e:
            logger.warning(
                "[INGEST][LITERATURE] Error detecting/ingesting signatures for item %s: %r",
                item_key,
                e,
            )
            # Non-blocking - continue

        logger.info("[INGEST][ZOTERO] Ingestion complete for item %s", item_key)
    else:
        logger.info(
            "[INGEST][ZOTERO] No new or changed attachments/notes; nothing ingested for item %s",
            item_key,
        )

