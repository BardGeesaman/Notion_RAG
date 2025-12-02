# amprenta_rag/ingestion/email_ingestion.py

"""
Email & Note ingestion into the Amprenta RAG engine.

Refactored from backup_2025_11_30/email_ingestion.py to:
- Use amprenta_rag.config + clients (OpenAI, Pinecone, Notion)
- Reuse the existing RAG DB schema:
    - Chunk ID (title)
    - Chunk Text (rich_text)
    - Parent Type (select)
    - Parent Email/Note Item (relation) for Email/Note
    - Parent Item (relation) for other types
    - Order (number)
    - Embedding Vector ID (rich_text)
    - Last Embedded (date)
- Operate on Notion Email DB where:
    - "Embedding Status" = "Not Embedded"
"""

from __future__ import annotations

import time
from datetime import datetime, timezone
from typing import Dict, Any, List, Optional

import requests

from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.zotero_ingest import _chunk_text, _embed_texts
from amprenta_rag.ingestion.notion_pages import (
    create_rag_chunk_page,
    fetch_not_embedded_emails,
    extract_page_content,
    update_email_page,
)
from amprenta_rag.ingestion.metadata_semantic import get_email_semantic_metadata

logger = get_logger(__name__)

REQUESTS_PER_SECOND = 2.0


# ----------------- Config helpers ----------------- #


def _notion_base_url() -> str:
    """Get Notion API base URL from config."""
    return get_config().notion.base_url


def _rag_db_id() -> str:
    """Get RAG DB ID from config."""
    cfg = get_config().notion
    if not cfg.rag_db_id:
        raise RuntimeError("Notion RAG DB ID is not configured in config.py")
    return cfg.rag_db_id

# ----------------- Cleanup orphaned chunks ----------------- #


def cleanup_orphaned_chunks() -> None:
    """
    Find and delete chunks whose parent email/note no longer exists.
    """
    logger.info("[INGEST][EMAIL] Checking for orphaned chunks...")

    db_id = _rag_db_id()
    query_url = f"{_notion_base_url()}/databases/{db_id}/query"

    payload: Dict[str, Any] = {
        "page_size": 100,
        "filter": {
            "or": [
                {
                    "property": "Parent Type",
                    "select": {"equals": "Email"},
                },
                {
                    "property": "Parent Type",
                    "select": {"equals": "Note"},
                },
            ]
        },
    }

    orphaned_chunks: List[Dict[str, Any]] = []
    chunk_ids_for_pinecone: List[str] = []
    next_cursor: Optional[str] = None

    while True:
        body = payload.copy()
        if next_cursor:
            body["start_cursor"] = next_cursor

        try:
            resp = requests.post(query_url, headers=notion_headers(), json=body)
            resp.raise_for_status()
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error querying RAG chunks: %r",
                e,
            )
            raise
        data = resp.json()

        results = data.get("results", [])
        if not results and not data.get("has_more"):
            break

        # Check each chunk's parent existence
        for chunk in results:
            props = chunk.get("properties", {})
            rel = props.get("Parent Email/Note Item", {}).get("relation", [])
            if not rel:
                # No parent at all → orphaned
                orphaned_chunks.append(chunk)
                continue

            parent_id = rel[0].get("id")
            if not parent_id:
                orphaned_chunks.append(chunk)
                continue

            check_url = f"{_notion_base_url()}/pages/{parent_id}"
            try:
                check_resp = requests.get(check_url, headers=notion_headers())
                if check_resp.status_code == 404:
                    orphaned_chunks.append(chunk)
            except Exception as e:
                logger.warning(
                    "[INGEST][EMAIL] Notion API error checking parent page %s: %r",
                    parent_id,
                    e,
                )
                # Treat as orphaned if we can't check
                orphaned_chunks.append(chunk)

            time.sleep(0.5 / REQUESTS_PER_SECOND)

        if not data.get("has_more"):
            break
        next_cursor = data.get("next_cursor")

        time.sleep(1.0 / REQUESTS_PER_SECOND)

    if not orphaned_chunks:
        logger.info("[INGEST][EMAIL] No orphaned chunks found")
        return

    logger.info("[INGEST][EMAIL] Deleting %d orphaned chunks...", len(orphaned_chunks))

    # Collect chunk page IDs + Pinecone IDs
    chunk_page_ids: List[str] = []
    for chunk in orphaned_chunks:
        cid_prop = chunk.get("properties", {}).get("Chunk ID", {}).get("title", [])
        if cid_prop:
            chunk_ids_for_pinecone.append(cid_prop[0].get("plain_text", ""))

        chunk_page_ids.append(chunk["id"].replace("-", ""))

    # Delete chunk pages from Notion
    for idx, chunk_page_id in enumerate(chunk_page_ids, 1):
        delete_url = f"{_notion_base_url()}/blocks/{chunk_page_id}"
        try:
            del_resp = requests.delete(delete_url, headers=notion_headers())
            if del_resp.status_code >= 300:
                logger.warning(
                    "[INGEST][EMAIL] Failed to delete chunk %d: %s",
                    idx,
                    del_resp.text,
                )
            elif idx % 10 == 0 or idx == len(chunk_page_ids):
                logger.info(
                    "[INGEST][EMAIL] Deleted %d/%d chunks",
                    idx,
                    len(chunk_page_ids),
                )
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error deleting chunk %s: %r",
                chunk_page_id,
                e,
            )
        time.sleep(1.0 / REQUESTS_PER_SECOND)

    # Delete vectors from Pinecone
    if chunk_ids_for_pinecone:
        logger.info(
            "[INGEST][EMAIL] Deleting %d vectors from Pinecone...",
            len(chunk_ids_for_pinecone),
        )
        index = get_pinecone_index()
        cfg = get_config()
        try:
            index.delete(ids=chunk_ids_for_pinecone, namespace=cfg.pinecone.namespace)
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Pinecone API error deleting orphaned vectors: %r",
                e,
            )
            raise

    logger.info("[INGEST][EMAIL] Cleaned up %d orphaned chunks", len(orphaned_chunks))


# ----------------- Delete email + chunks ----------------- #


def delete_email_and_chunks(email_page_id: str) -> None:
    """
    Delete an email and all its associated RAG chunks from Notion and Pinecone.
    """
    try:
        logger.info("[INGEST][EMAIL] Starting deletion for email: %s", email_page_id)

        # 1. Query RAG Engine for all chunks with this email as parent
        db_id = _rag_db_id()
        query_url = f"{_notion_base_url()}/databases/{db_id}/query"

        payload = {
            "filter": {
                "property": "Parent Email/Note Item",
                "relation": {"contains": email_page_id},
            }
        }

        try:
            resp = requests.post(query_url, headers=notion_headers(), json=payload)
            if resp.status_code >= 300:
                logger.error("[INGEST][EMAIL] Error querying RAG chunks: %s", resp.text)
                resp.raise_for_status()
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error querying RAG chunks for email %s: %r",
                email_page_id,
                e,
            )
            raise

        chunks = resp.json().get("results", [])
        logger.info("[INGEST][EMAIL] Found %d chunks to delete for email %s", len(chunks), email_page_id)

        chunk_ids: List[str] = []
        chunk_page_ids: List[str] = []
        for chunk in chunks:
            props = chunk.get("properties", {})
            cid_prop = props.get("Chunk ID", {}).get("title", [])
            if cid_prop:
                chunk_ids.append(cid_prop[0].get("plain_text", ""))
            chunk_page_ids.append(chunk["id"].replace("-", ""))

        # 2. Delete chunk pages from Notion
        for idx, chunk_page_id in enumerate(chunk_page_ids, 1):
            delete_url = f"{_notion_base_url()}/blocks/{chunk_page_id}"
            try:
                del_resp = requests.delete(delete_url, headers=notion_headers())
                if del_resp.status_code >= 300:
                    logger.warning(
                        "[INGEST][EMAIL] Failed to delete chunk page %s: %s",
                        chunk_page_id,
                        del_resp.text,
                    )
                elif idx % 10 == 0 or idx == len(chunk_page_ids):
                    logger.info(
                        "[INGEST][EMAIL] Deleted %d/%d chunks",
                        idx,
                        len(chunk_page_ids),
                    )
            except Exception as e:
                logger.error(
                    "[INGEST][EMAIL] Notion API error deleting chunk page %s: %r",
                    chunk_page_id,
                    e,
                )

            time.sleep(1.0 / REQUESTS_PER_SECOND)

        # 3. Delete vectors from Pinecone
        if chunk_ids:
            logger.info(
                "[INGEST][EMAIL] Deleting %d vectors from Pinecone for email %s...",
                len(chunk_ids),
                email_page_id,
            )
            index = get_pinecone_index()
            cfg = get_config()
            try:
                index.delete(ids=chunk_ids, namespace=cfg.pinecone.namespace)
                logger.info("[INGEST][EMAIL] Pinecone vectors deleted for email %s", email_page_id)
            except Exception as e:
                logger.error(
                    "[INGEST][EMAIL] Pinecone API error deleting vectors for email %s: %r",
                    email_page_id,
                    e,
                )
                raise
        else:
            logger.info("[INGEST][EMAIL] No vectors to delete from Pinecone for email %s", email_page_id)

        # 4. Delete the email page itself
        logger.info("[INGEST][EMAIL] Deleting email page %s...", email_page_id)
        delete_url = f"{_notion_base_url()}/blocks/{email_page_id}"
        try:
            email_resp = requests.delete(delete_url, headers=notion_headers())
            if email_resp.status_code >= 300:
                logger.error(
                    "[INGEST][EMAIL] Failed to delete email page %s: %s",
                    email_page_id,
                    email_resp.text,
                )
                raise Exception(f"Failed to delete email page: {email_resp.text}")
            else:
                logger.info("[INGEST][EMAIL] Email page deleted: %s", email_page_id)
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error deleting email page %s: %r",
                email_page_id,
                e,
            )
            raise

        logger.info("[INGEST][EMAIL] Deletion complete for email %s", email_page_id)

    except Exception as e:
        logger.error("[INGEST][EMAIL] Deletion failed for email %s: %r", email_page_id, e)
        raise


# ----------------- Main ingestion ----------------- #


def ingest_email(email_page: Dict[str, Any], parent_type: str = "Email") -> None:
    """
    Ingest a single email/note page into RAG Engine + Pinecone.
    
    This function:
    1. Extracts metadata (title, from, tags, semantic metadata) from the email page
    2. Extracts full text content from Notion blocks
    3. Chunks the text using shared chunking logic
    4. Embeds chunks using OpenAI
    5. Creates RAG chunk pages in Notion
    6. Upserts vectors to Pinecone with rich metadata
    7. Updates the email page status to "Embedded"
    
    Metadata structure matches zotero_ingest.py for consistency:
    - Doc-level fields: doc_id, doc_source, doc_type, diseases, targets, etc.
    - Lipid-level fields: lipids, lipid_classes, lipid_signatures, etc.
    - Chunk-level fields: chunk_id, chunk_index, snippet, etc.
    
    Args:
        email_page: Notion page object from Email DB query
        parent_type: Type label ("Email" or "Note") for the parent
    """
    index = get_pinecone_index()
    cfg = get_config()

    page_id = email_page["id"].replace("-", "")
    props = email_page.get("properties", {})

    # Extract metadata
    title_prop = props.get("Title", {}).get("title", [])
    title = title_prop[0]["plain_text"] if title_prop else "Untitled"

    from_prop = props.get("From", {}).get("rich_text", [])
    from_sender = from_prop[0]["plain_text"] if from_prop else ""

    # Extract tags from multi-select property "Tags" (if present)
    tags_prop = props.get("Tags", {})
    tags: List[str] = []
    if "multi_select" in tags_prop:
        tags = [t["name"] for t in tags_prop["multi_select"] if t.get("name")]

    type_prop = props.get("Type", {}).get("select", {})
    item_type = type_prop.get("name", "Email")

    # Doc-level + lipid-level metadata from Email DB
    semantic_meta = get_email_semantic_metadata(email_page)

    # Identity / doc-level fields
    doc_meta: Dict[str, Any] = {
        "doc_id": f"EMAIL:{page_id}",
        "doc_source": "Email",
        "doc_source_subtype": "EmailBody",
        "doc_type": "InternalNote" if item_type == "Note" else "Email",
        "diseases": semantic_meta.get("diseases", []),
        "targets": semantic_meta.get("targets", []),
        "modality": semantic_meta.get("modality", []),
        "stage": semantic_meta.get("stage"),
        "model_systems": semantic_meta.get("model_systems", []),
        "biomarker_role": semantic_meta.get("biomarker_role", []),
        "importance": semantic_meta.get("importance"),
        "manual_tags": tags,
        # lipid-level
        "lipids_raw": semantic_meta.get("lipids_raw", []),
        "lipids": semantic_meta.get("lipids", []),
        "lipid_classes": semantic_meta.get("lipid_classes", []),
        "lipid_signatures": semantic_meta.get("lipid_signatures", []),
        "lipid_signature_role": semantic_meta.get("lipid_signature_role", []),
        "phenotype_axes": semantic_meta.get("phenotype_axes", []),
        "signature_ownership": semantic_meta.get("signature_ownership", []),
        "matrix": semantic_meta.get("matrix", []),
        "treatment_arms": semantic_meta.get("treatment_arms", []),
    }    

    logger.info("[INGEST][EMAIL] Processing: %s (%s)", title[:60], item_type)

    # Extract full content from blocks
    try:
        full_text = extract_page_content(page_id)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Error extracting page content for email %s: %r",
            page_id,
            e,
        )
        raise

    # Add metadata header for more context
    header = f"Title: {title}\n"
    if from_sender:
        header += f"From: {from_sender}\n"
    header += f"Type: {item_type}\n\n"

    full_text = header + full_text

    if len(full_text.strip()) < 50:
        logger.info(
            "[INGEST][EMAIL] Skipping - content too short (%d chars) for email %s",
            len(full_text),
            page_id,
        )
        return

    logger.info("[INGEST][EMAIL] Extracted %d characters for email %s", len(full_text), page_id)

    # Chunk and embed using shared helpers from zotero_ingest
    chunks = _chunk_text(full_text)
    logger.info("[INGEST][EMAIL] Generated %d chunks for email %s", len(chunks), page_id)

    if not chunks:
        logger.info("[INGEST][EMAIL] No chunks to embed for email %s", page_id)
        return

    logger.info("[INGEST][EMAIL] Embedding chunks with OpenAI for email %s", page_id)
    try:
        embeddings = _embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] OpenAI API error embedding chunks for email %s: %r",
            page_id,
            e,
        )
        raise

    now = datetime.now(timezone.utc).isoformat()
    vectors = []

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id}_chunk_{order:03d}"

        # Create Notion chunk page using shared helper from notion_pages
        # The helper automatically handles Email/Note relation property selection
        try:
            chunk_page_id = create_rag_chunk_page(
                chunk_id, chunk, parent_type, page_id, order, when_iso=now
            )
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error creating chunk page %s: %r",
                chunk_id,
                e,
            )
            raise

        snippet = chunk[:200]

        # Build metadata for Pinecone
        metadata: Dict[str, Any] = {
            **doc_meta,      # doc-level + lipid-level from Email DB
            "chunk_id": chunk_id,
            "chunk_index": order,
            "section": "Body",
            "snippet": snippet,
            "notion_chunk_page_id": chunk_page_id,
            "parent_email_page_id": page_id,
            "source": "Email",
            "source_type": parent_type,  # "Email" or "Note"
            "title": title[:200],
            "from": (from_sender[:200] if from_sender else ""),
            "tags": tags,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(metadata),
            }
        )

    logger.info(
        "[INGEST][EMAIL] Upserting %d vectors into Pinecone for email %s...",
        len(vectors),
        page_id,
    )
    try:
        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Pinecone API error upserting vectors for email %s: %r",
            page_id,
            e,
        )
        raise

    logger.info("[INGEST][EMAIL] Updating Email page status for email %s...", page_id)
    try:
        update_email_page(page_id, now)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Notion API error updating email page %s: %r",
            page_id,
            e,
        )
        raise

    logger.info("[INGEST][EMAIL] Email ingestion complete for email %s", page_id)


def batch_ingest_emails(parent_type: str = "Email") -> None:
    """
    Batch process all not-embedded emails from the Email DB.
    
    This function:
    1. Cleans up orphaned chunks (chunks whose parent emails no longer exist)
    2. Fetches all emails with "Embedding Status" = "Not Embedded"
    3. Processes each email through ingest_email()
    4. Handles errors gracefully, continuing with remaining emails
    
    Args:
        parent_type: Type label to use for ingested emails ("Email" or "Note")
    """
    cleanup_orphaned_chunks()

    logger.info("[INGEST][EMAIL] " + "=" * 60)

    try:
        emails = fetch_not_embedded_emails()
    except Exception as e:
        logger.error("[INGEST][EMAIL] Error fetching not-embedded emails: %r", e)
        raise

    if not emails:
        logger.info("[INGEST][EMAIL] All emails are already embedded!")
        return

    total = len(emails)
    logger.info("[INGEST][EMAIL] Starting batch ingestion of %d email(s)...", total)

    for idx, email_page in enumerate(emails, 1):
        props = email_page.get("properties", {})
        title_prop = props.get("Title", {}).get("title", [])
        title = title_prop[0]["plain_text"] if title_prop else "Untitled"

        logger.info("[INGEST][EMAIL] [%d/%d] %s", idx, total, title[:80])

        try:
            ingest_email(email_page, parent_type=parent_type)
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Error processing email [%d/%d] %s: %r",
                idx,
                total,
                title[:80],
                e,
            )

        time.sleep(1.0 / REQUESTS_PER_SECOND)

    logger.info("[INGEST][EMAIL] Batch complete: processed %d email(s)", total)