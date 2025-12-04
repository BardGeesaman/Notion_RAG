"""
Attachment processing for Zotero ingestion.

This module provides functions for processing Zotero attachments:
downloading, extracting text, chunking, embedding, and creating Notion pages.
"""

from __future__ import annotations

import textwrap
from datetime import datetime, timezone
from typing import Any, Dict, List

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.notion_pages import create_rag_chunk_page
from amprenta_rag.ingestion.pinecone_utils import (
    attachment_already_ingested,
    sanitize_metadata,
)
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.ingestion.text_extraction import (
    extract_text_from_attachment_bytes,
    is_boilerplate,
)
from amprenta_rag.ingestion.zotero_api import (
    ZoteroItem,
    download_zotero_file,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def process_attachments(
    item: ZoteroItem,
    attachments: List[Dict[str, Any]],
    parent_page_id: str,
    parent_type: str,
    base_meta: Dict[str, Any],
    force: bool = False,
) -> tuple[List[str], bool]:
    """
    Process all attachments for a Zotero item.

    Args:
        item: Zotero item metadata
        attachments: List of attachment dictionaries from Zotero API
        parent_page_id: Notion page ID of the parent literature item
        parent_type: Type of parent page (e.g., "Literature")
        base_meta: Base metadata dictionary to include in all vectors
        force: If True, re-ingest even if already ingested

    Returns:
        Tuple of (list of extracted text strings, whether any were ingested)
    """
    cfg = get_config()
    index = get_pinecone_index()
    now = datetime.now(timezone.utc).isoformat()
    any_ingested = False
    all_text_parts: List[str] = []

    for att in attachments:
        att_key = att.get("key")
        if not att_key:
            continue

        filename = att.get("filename") or ""
        content_type = att.get("contentType")
        att_md5 = att.get("md5")

        logger.info(
            "[INGEST][ZOTERO] Processing attachment %s (%s)",
            att_key,
            filename or "unnamed",
        )

        # Skip Zotero HTML snapshots
        title = att.get("title", "") or ""
        if "snapshot" in filename.lower() or "snapshot" in title.lower():
            logger.info("[INGEST][ZOTERO] Skipping Snapshot attachment %s", att_key)
            continue
        if content_type and content_type.lower().startswith("text/html"):
            logger.info(
                "[INGEST][ZOTERO] Skipping HTML snapshot attachment %s", att_key
            )
            continue

        if not force and attachment_already_ingested(
            item_key=item.key,
            attachment_key=att_key,
            current_md5=att_md5,
        ):
            logger.info(
                "[INGEST][ZOTERO] Already ingested with matching md5; skipping attachment %s",
                att_key,
            )
            continue

        try:
            data = download_zotero_file(att_key)
        except Exception as e:
            logger.error(
                "[INGEST][ZOTERO] Zotero API error downloading attachment %s: %r",
                att_key,
                e,
            )
            raise

        text = extract_text_from_attachment_bytes(content_type, data, filename=filename)
        if not text or len(text.strip()) < 100:
            logger.info(
                "[INGEST][ZOTERO] Attachment text is empty or very short; skipping %s",
                att_key,
            )
            continue

        header = item.title + "\n"
        if item.journal:
            header += f"{item.journal}\n"
        if item.doi:
            header += f"DOI: {item.doi}\n"
        header += f"[Attachment: {filename}]\n\n"
        full_text = header + text

        # Collect text for feature extraction
        all_text_parts.append(text)

        chunks = chunk_text(full_text)
        chunks = [c for c in chunks if not is_boilerplate(c)]
        if not chunks:
            logger.info(
                "[INGEST][ZOTERO] No usable chunks produced; skipping attachment %s",
                att_key,
            )
            continue

        logger.info(
            "[INGEST][ZOTERO] Generated %d chunk(s) for attachment %s",
            len(chunks),
            att_key,
        )
        logger.info(
            "[INGEST][ZOTERO] Embedding chunks with OpenAI for attachment %s", att_key
        )
        embeddings = embed_texts(chunks)

        vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{parent_page_id}_{att_key}_chunk_{order:03d}"

            try:
                chunk_page_id = create_rag_chunk_page(
                    chunk_id=chunk_id,
                    chunk_text=chunk,
                    parent_type=parent_type,
                    parent_page_id=parent_page_id,
                    order=order,
                    when_iso=now,
                )
            except Exception as e:
                logger.error(
                    "[INGEST][ZOTERO] Notion API error creating chunk page %s: %r",
                    chunk_id,
                    e,
                )
                raise

            snippet = textwrap.shorten(chunk, width=300)

            meta: Dict[str, Any] = {
                **base_meta,
                "chunk_id": chunk_id,
                "chunk_index": order,
                "snippet": snippet,
                "zotero_item_key": item.key,
                "attachment_key": att_key,
                "attachment_filename": filename,
                "attachment_content_type": content_type,
                "attachment_md5": att_md5,
                "notion_chunk_page_id": chunk_page_id,
                "source": "Literature",
                "source_type": parent_type,
                "title": item.title,
                "zotero_tags": item.tags,
                "item_type": item.item_type,
            }
            if item.doi:
                meta["doi"] = item.doi
            if item.url:
                meta["url"] = item.url
            if item.date:
                meta["date"] = item.date

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        if vectors:
            logger.info(
                "[INGEST][ZOTERO] Upserting %d vectors into Pinecone for attachment %s",
                len(vectors),
                att_key,
            )
            try:
                index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
            except Exception as e:
                logger.error(
                    "[INGEST][ZOTERO] Pinecone API error upserting vectors for attachment %s: %r",
                    att_key,
                    e,
                )
                raise
            any_ingested = True

    return all_text_parts, any_ingested

