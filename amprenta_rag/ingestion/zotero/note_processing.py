"""
Note processing for Zotero ingestion.

This module provides functions for processing Zotero notes:
extracting text, chunking, embedding, and creating Notion pages.
"""

from __future__ import annotations

import hashlib
import textwrap
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.ingestion.notion_pages import create_rag_chunk_page
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

def create_rag_chunk_page(chunk_id: str, chunk_text: str, parent_type: str, parent_id: str, order: int, when_iso: str) -> Optional[str]:
    """DEPRECATED: Notion support removed. Returns None."""
    logger.debug("[ZOTERO][NOTE] create_rag_chunk_page() deprecated - Notion support removed")
    return None
from amprenta_rag.ingestion.pinecone_utils import note_already_ingested, sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.ingestion.text_extraction import html_to_text, is_boilerplate
from amprenta_rag.ingestion.zotero_api import ZoteroItem
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def process_notes(
    item: ZoteroItem,
    notes: List[Dict[str, Any]],
    parent_page_id: str,
    parent_type: str,
    base_meta: Dict[str, Any],
    force: bool = False,
) -> tuple[List[str], bool]:
    """
    Process all notes for a Zotero item.

    Args:
        item: Zotero item metadata
        notes: List of note dictionaries from Zotero API
        parent_page_id: Notion page ID of the parent literature item
        parent_type: Type of parent page (e.g., "Literature")
        base_meta: Base metadata dictionary to include in all vectors
        force: If True, re-ingest even if already ingested

    Returns:
        Tuple of (list of extracted text strings, whether any were ingested)
    """
    cfg = get_config()
    store = get_vector_store()
    now = datetime.now(timezone.utc).isoformat()
    any_ingested = False
    all_text_parts: List[str] = []

    for note in notes:
        note_key = note.get("key")
        if not note_key:
            continue

        raw_html = note.get("note") or ""
        text = html_to_text(raw_html)
        if not text or len(text.strip()) < 50:
            continue

        note_hash = hashlib.md5(text.encode("utf-8", errors="ignore")).hexdigest()
        logger.info(
            "[INGEST][ZOTERO] Processing note %s (hash=%s)", note_key, note_hash
        )

        # Collect text for feature extraction
        all_text_parts.append(text)

        if not force and note_already_ingested(
            item_key=item.key,
            note_key=note_key,
            note_hash=note_hash,
        ):
            logger.info(
                "[INGEST][ZOTERO] Note already ingested with same hash; skipping %s",
                note_key,
            )
            continue

        header = f"{item.title}\n[Zotero Note]\n\n"
        full_text = header + text

        chunks = chunk_text(full_text)
        chunks = [c for c in chunks if not is_boilerplate(c)]
        if not chunks:
            logger.info(
                "[INGEST][ZOTERO] No usable chunks for note; skipping %s", note_key
            )
            continue

        logger.info(
            "[INGEST][ZOTERO] Generated %d chunk(s) for note %s",
            len(chunks),
            note_key,
        )
        logger.info(
            "[INGEST][ZOTERO] Embedding note chunks with OpenAI for note %s", note_key
        )
        embeddings = embed_texts(chunks)

        note_vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{parent_page_id}_{note_key}_note_{order:03d}"

            try:
                chunk_page_id = create_rag_chunk_page(
                    chunk_id=chunk_id,
                    chunk_text=chunk,
                    parent_type=parent_type,
                    parent_id=parent_page_id,
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

            note_meta: Dict[str, Any] = {
                **base_meta,
                "chunk_id": chunk_id,
                "chunk_index": order,
                "snippet": snippet,
                "zotero_item_key": item.key,
                "note_key": note_key,
                "note_hash": note_hash,
                "notion_chunk_page_id": chunk_page_id,
                "source": "Literature",
                "source_type": parent_type,
                "title": item.title,
                "zotero_tags": item.tags,
                "item_type": item.item_type,
            }
            if item.doi:
                note_meta["doi"] = item.doi
            if item.url:
                note_meta["url"] = item.url
            if item.date:
                note_meta["date"] = item.date

            note_vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(note_meta),
                }
            )

        if note_vectors:
            logger.info(
                "[INGEST][ZOTERO] Upserting %d vectors into Pinecone for note %s",
                len(note_vectors),
                note_key,
            )
            try:
                store.upsert(vectors=note_vectors, namespace=cfg.pinecone.namespace)
            except Exception as e:
                logger.error(
                    "[INGEST][ZOTERO] Pinecone API error upserting vectors for note %s: %r",
                    note_key,
                    e,
                )
                raise
            any_ingested = True

    return all_text_parts, any_ingested

