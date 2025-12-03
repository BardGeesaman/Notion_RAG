"""
Email cleanup utilities.

This module handles cleanup operations for email/note ingestion:
- Finding and deleting orphaned chunks
- Deleting emails and their associated chunks
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

REQUESTS_PER_SECOND = 2.0

__all__ = [
    "cleanup_orphaned_chunks",
    "delete_email_and_chunks",
]


def _notion_base_url() -> str:
    """Get Notion API base URL from config."""
    return get_config().notion.base_url


def _rag_db_id() -> str:
    """Get RAG DB ID from config."""
    cfg = get_config().notion
    if not cfg.rag_db_id:
        raise RuntimeError("Notion RAG DB ID is not configured in config.py")
    return cfg.rag_db_id


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
            resp = requests.post(
                query_url, headers=notion_headers(), json=body, timeout=30
            )
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
                check_resp = requests.get(
                    check_url, headers=notion_headers(), timeout=30
                )
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
            del_resp = requests.delete(
                delete_url, headers=notion_headers(), timeout=30
            )
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
            resp = requests.post(
                query_url, headers=notion_headers(), json=payload, timeout=30
            )
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
        logger.info(
            "[INGEST][EMAIL] Found %d chunks to delete for email %s",
            len(chunks),
            email_page_id,
        )

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
                del_resp = requests.delete(
                    delete_url, headers=notion_headers(), timeout=30
                )
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
                logger.info(
                    "[INGEST][EMAIL] Pinecone vectors deleted for email %s",
                    email_page_id,
                )
            except Exception as e:
                logger.error(
                    "[INGEST][EMAIL] Pinecone API error deleting vectors for email %s: %r",
                    email_page_id,
                    e,
                )
                raise
        else:
            logger.info(
                "[INGEST][EMAIL] No vectors to delete from Pinecone for email %s",
                email_page_id,
            )

        # 4. Delete the email page itself
        logger.info("[INGEST][EMAIL] Deleting email page %s...", email_page_id)
        delete_url = f"{_notion_base_url()}/blocks/{email_page_id}"
        try:
            email_resp = requests.delete(
                delete_url, headers=notion_headers(), timeout=30
            )
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
        logger.error(
            "[INGEST][EMAIL] Deletion failed for email %s: %r", email_page_id, e
        )
        raise

