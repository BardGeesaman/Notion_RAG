# amprenta_rag/ingestion/zotero_collection.py

from __future__ import annotations

from typing import Dict, Any, List, Iterable, Set, Optional

import requests

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index


from amprenta_rag.ingestion.zotero_ingest import ingest_zotero_item

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
    resp = requests.get(url, headers=_zotero_headers())
    resp.raise_for_status()
    items = resp.json()
    if not isinstance(items, list):
        raise RuntimeError(f"Unexpected Zotero response: {items!r}")
    return items


# ---------- Notion helpers ---------- #


def _require_notion_db_id(name: str, value: Optional[str]) -> str:
    if not value:
        raise RuntimeError(
            f"Notion {name} DB ID is not configured. "
            f"Set AMPRENTA_NOTION_{name}_DB_ID."
        )
    return value


def get_literature_page_id_for_item(item_key: str) -> Optional[str]:
    """
    Find the Literature page whose 'Zotero Item Key' == item_key.
    """
    cfg = get_config().notion
    lit_db_id = _require_notion_db_id("LIT", cfg.lit_db_id)

    url = f"{cfg.base_url}/databases/{lit_db_id}/query"
    payload = {
        "page_size": 1,
        "filter": {
            "property": "Zotero Item Key",
            "rich_text": {"equals": item_key},
        },
    }

    resp = requests.post(url, headers=notion_headers(), json=payload)
    if resp.status_code >= 300:
        logger.warning(
            "Notion Literature query error for item %s: %s %s",
            item_key,
            resp.status_code,
            resp.text,
        )
        resp.raise_for_status()

    results = resp.json().get("results", [])
    if not results:
        return None

    return results[0]["id"]


def delete_rag_pages_for_item(item_key: str) -> None:
    """
    Archive (soft-delete) all RAG chunk pages whose "Parent Item" points
    to the Literature page for this Zotero item.
    """
    cfg = get_config().notion
    rag_db_id = _require_notion_db_id("RAG", cfg.rag_db_id)

    lit_page_id = get_literature_page_id_for_item(item_key)
    if lit_page_id is None:
        logger.info("   ℹ️ No Literature page; skipping RAG page deletion.")
        return

    url = f"{cfg.base_url}/databases/{rag_db_id}/query"
    payload = {
        "page_size": 100,
        "filter": {
            "property": "Parent Item",
            "relation": {"contains": lit_page_id},
        },
    }

    resp = requests.post(url, headers=notion_headers(), json=payload)
    if resp.status_code >= 300:
        logger.warning(
            "Notion RAG query error for item %s: %s %s",
            item_key,
            resp.status_code,
            resp.text,
        )
        resp.raise_for_status()

    results = resp.json().get("results", [])
    if not results:
        logger.info("   ℹ️ No RAG pages found.")
        return

    for page in results:
        page_id = page["id"]
        logger.info("   🗑️ Archiving RAG page %s", page_id.replace("-", ""))
        requests.patch(
            f"{cfg.base_url}/pages/{page_id}",
            headers=notion_headers(),
            json={"archived": True},
        )


# ---------- Pinecone helpers ---------- #


def delete_pinecone_vectors_for_item(item_key: str) -> None:
    """
    Delete all Pinecone vectors whose metadata.zotero_item_key == item_key.
    """
    index = get_pinecone_index()
    index.delete(filter={"zotero_item_key": {"$eq": item_key}})


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
            logger.info(
                "⏭️ Skipping non-literature item %s [%s]", item_key, item_type
            )
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

        print(f"\n=== [{idx}/{total}] {title} ({item_key}) [{item_type}] ===")

        try:
            ingest_zotero_item(item_key=item_key, parent_type=parent_type)
        except Exception as e:  # noqa: BLE001
            print(f"❌ Skipping item {item_key} due to error: {e}")
            continue

        print("")

    print("\n🎉 Incremental collection ingest complete.")


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

        print(f"\n=== [{idx}/{total}] {title} ({item_key}) [{item_type}] ===")

        print(f"🧹 Clearing old data for {item_key}...")
        delete_pinecone_vectors_for_item(item_key)
        delete_rag_pages_for_item(item_key)

        print(f"🔄 Re-ingesting {item_key}...")
        try:
            ingest_zotero_item(
                item_key=item_key,
                parent_type=parent_type,
                force=force,
            )

        except Exception as e:  # noqa: BLE001
            print(f"❌ Skipping item {item_key} due to error: {e}")
            continue

        print("")

    print("\n🎉 Collection resync complete.")