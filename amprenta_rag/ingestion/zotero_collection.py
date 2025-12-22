# amprenta_rag/ingestion/zotero_collection.py
"""
Zotero collection ingestion utilities.

Handles ingestion of Zotero collections into Pinecone.
Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set

import requests

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.zotero_ingest import ingest_zotero_item
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Default set of "literature-like" item types
DEFAULT_ALLOWED_TYPES: Set[str] = {
    "journalArticle",
    "book",
    "bookSection",
    "report",
    "preprint",
    "conferencePaper",
}


# ---------- Zotero helpers ---------- #


def _zotero_headers() -> Dict[str, str]:
    cfg = get_config().zotero
    if not cfg.api_key:
        raise RuntimeError("AMPRENTA_ZOTERO_API_KEY is not set.")
    return {
        "Zotero-API-Key": cfg.api_key,
        "Accept": "application/json",
    }


def _zotero_base() -> str:
    cfg = get_config().zotero
    if not cfg.library_type or not cfg.library_id:
        raise RuntimeError(
            "Zotero library not configured. "
            "Set AMPRENTA_ZOTERO_LIBRARY_TYPE and AMPRENTA_ZOTERO_LIBRARY_ID."
        )
    return f"https://api.zotero.org/{cfg.library_type}s/{cfg.library_id}"


def fetch_collection_items(collection_key: str) -> List[Dict[str, Any]]:
    """
    Fetch top-level (parent) items for a given collection.

    Mirrors the behavior of your previous ingest_collection.py and resync_collection.py,
    using the /items/top endpoint with include=data.
    """
    url = f"{_zotero_base()}/collections/{collection_key}/items/top?include=data"
    resp = requests.get(url, headers=_zotero_headers(), timeout=30)
    resp.raise_for_status()
    items = resp.json()
    if not isinstance(items, list):
        raise RuntimeError(f"Unexpected Zotero response: {items!r}")
    return items


# ---------- Notion helpers (stubs - Notion support removed) ---------- #


def get_literature_page_id_for_item(item_key: str) -> Optional[str]:
    """
    Stub: Notion support removed. Returns None.

    Previously found Literature page by Zotero Item Key in Notion.
    """
    logger.debug("[ZOTERO-COLLECTION] get_literature_page_id_for_item() is a no-op (Notion removed)")
    return None


def delete_rag_pages_for_item(item_key: str) -> None:
    """
    Stub: Notion support removed. No-op.

    Previously archived RAG chunk pages in Notion.
    """
    logger.debug("[ZOTERO-COLLECTION] delete_rag_pages_for_item() is a no-op (Notion removed)")


# ---------- Pinecone helpers ---------- #


def delete_pinecone_vectors_for_item(item_key: str) -> None:
    """
    Delete all Pinecone vectors whose metadata.zotero_item_key == item_key.
    """
    cfg = get_config()
    store = get_vector_store()
    store.delete(filter={"zotero_item_key": {"$eq": item_key}}, namespace=cfg.pinecone.namespace)


# ---------- High-level collection operations ---------- #


def _iter_eligible_items(
    items: Iterable[Dict[str, Any]],
    allowed_types: Optional[Set[str]] = None,
) -> Iterable[Dict[str, Any]]:
    if allowed_types is None:
        allowed_types = DEFAULT_ALLOWED_TYPES

    for item in items:
        data = item.get("data", {}) or {}
        item_key = data.get("key")
        item_type = data.get("itemType")

        if not item_key:
            continue

        if item_type not in allowed_types:
            logger.info("⏭️ Skipping non-literature item %s [%s]", item_key, item_type)
            continue

        yield item


def incremental_ingest_collection(
    collection_key: str,
    *,
    parent_type: str = "Literature",
    allowed_types: Optional[Set[str]] = None,
) -> None:
    """
    Incrementally ingest all *eligible* items in a Zotero collection.

    For each eligible parent item:
      - Call ingest_zotero_item(item_key, parent_type)

    Because ingest_zotero_item is idempotent:
      - New items     -> fully ingested
      - Existing same -> skipped
      - New attachments/notes -> only new ones ingested
    """
    logger.info("📚 Fetching parent items for collection %s...", collection_key)
    items = fetch_collection_items(collection_key)
    total = len(items)
    logger.info("Found %d parent items.", total)

    for idx, item in enumerate(
        _iter_eligible_items(items, allowed_types=allowed_types), start=1
    ):
        data = item.get("data", {}) or {}
        item_key = data.get("key")
        title = data.get("title") or "(untitled)"
        item_type = data.get("itemType")

        logger.info("=== [%d/%d] %s (%s) [%s] ===", idx, total, title, item_key, item_type)

        try:
            if not isinstance(item_key, str):
                logger.warning("Skipping item with missing key in collection %s", collection_key)
                continue
            ingest_zotero_item(item_key=item_key, parent_type=parent_type)
        except Exception as e:  # noqa: BLE001
            logger.error("Skipping item %s due to error: %r", item_key, e)
            continue
        logger.info("")

    logger.info("Incremental collection ingest complete.")


def resync_collection(
    collection_key: str,
    *,
    parent_type: str = "Literature",
    allowed_types: Optional[Set[str]] = None,
    force: bool = False,
) -> None:
    """
    Fully resync a Zotero collection:

    For each eligible parent item:
      - Delete Pinecone vectors for this zotero_item_key
      - Delete RAG chunk pages in Notion (by Parent Item relation)
      - Re-ingest the item using ingest_zotero_item()
    """
    logger.info("📚 Fetching parent items for collection %s...", collection_key)
    items = fetch_collection_items(collection_key)
    total = len(items)
    logger.info("Found %d parent items.", total)

    for idx, item in enumerate(
        _iter_eligible_items(items, allowed_types=allowed_types), start=1
    ):
        data = item.get("data", {}) or {}
        item_key = data.get("key")
        title = data.get("title") or "(untitled)"
        item_type = data.get("itemType")

        logger.info("=== [%d/%d] %s (%s) [%s] ===", idx, total, title, item_key, item_type)

        if not isinstance(item_key, str):
            logger.warning("Skipping item with missing key in collection %s", collection_key)
            continue

        logger.info("🧹 Clearing old data for %s...", item_key)
        delete_pinecone_vectors_for_item(item_key)
        delete_rag_pages_for_item(item_key)

        logger.info("🔄 Re-ingesting %s...", item_key)
        try:
            ingest_zotero_item(
                item_key=item_key,
                parent_type=parent_type,
                force=force,
            )

        except Exception as e:  # noqa: BLE001
            logger.error("Skipping item %s due to error: %r", item_key, e)
            continue
        logger.info("")

    logger.info("Collection resync complete.")
