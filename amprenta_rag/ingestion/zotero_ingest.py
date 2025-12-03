# amprenta_rag/ingestion/zotero_ingest.py

"""
Zotero item ingestion pipeline.

This module orchestrates the complete ingestion of Zotero items:
1. Fetch Zotero item metadata and attachments/notes
2. Extract text from PDFs and notes
3. Chunk and embed text using OpenAI
4. Create Notion pages for literature items and chunks
5. Upsert vectors to Pinecone with rich metadata

The pipeline is idempotent - unchanged attachments/notes are skipped.
"""

from __future__ import annotations

import hashlib
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

from amprenta_rag.clients.openai_client import (get_default_models,
                                                get_openai_client)
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text, link_features_to_notion_items)
from amprenta_rag.ingestion.metadata_semantic import \
    get_literature_semantic_metadata
from amprenta_rag.ingestion.notion_pages import (create_rag_chunk_page,
                                                 ensure_literature_page,
                                                 update_literature_page)
from amprenta_rag.ingestion.pinecone_utils import (attachment_already_ingested,
                                                   note_already_ingested,
                                                   sanitize_metadata)
from amprenta_rag.ingestion.signature_integration import \
    detect_and_ingest_signatures_from_content
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.ingestion.text_extraction import (
    extract_text_from_attachment_bytes, html_to_text, is_boilerplate)
from amprenta_rag.ingestion.zotero_api import (ZoteroItem,
                                               download_zotero_file,
                                               fetch_zotero_attachments,
                                               fetch_zotero_item,
                                               fetch_zotero_notes)
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
    cfg = get_config()
    index = get_pinecone_index()

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
    any_ingested = False

    # Collect all text content for feature extraction
    all_text_parts: List[str] = []

    # ---------------- Attachments ---------------- #
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
            item_key=item_key,
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

    # ---------------- Notes ---------------- #
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
            item_key=item_key,
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
                meta["doi"] = item.doi
            if item.url:
                meta["url"] = item.url
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
                index.upsert(vectors=note_vectors, namespace=cfg.pinecone.namespace)
            except Exception as e:
                logger.error(
                    "[INGEST][ZOTERO] Pinecone API error upserting vectors for note %s: %r",
                    note_key,
                    e,
                )
                raise
            any_ingested = True

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
                from pathlib import Path

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
