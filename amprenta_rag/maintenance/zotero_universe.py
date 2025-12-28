# amprenta_rag/maintenance/zotero_universe.py
"""
Zotero universe maintenance utilities.

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from typing import Any, Dict, List, Set

import requests

from amprenta_rag.config import get_config
from amprenta_rag.ingestion import (incremental_ingest_collection,
                                    resync_collection)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


# --------- Zotero helpers --------- #


def _zotero_headers() -> Dict[str, str]:
    cfg = get_config().zotero
    return {
        "Zotero-API-Key": cfg.api_key,
        "Accept": "application/json",
    }


def _zotero_base() -> str:
    cfg = get_config().zotero
    return f"https://api.zotero.org/{cfg.library_type}s/{cfg.library_id}"


def zotero_item_exists(item_key: str) -> bool:
    """
    Return True if Zotero /items/{item_key} returns 200, False if 404.
    """
    url = f"{_zotero_base()}/items/{item_key}"
    resp = requests.get(url, headers=_zotero_headers())
    if resp.status_code == 404:
        return False
    resp.raise_for_status()
    return True


def collection_item_keys(collection_key: str) -> Set[str]:
    """
    Return the set of parent item keys in a Zotero collection.
    Mirrors legacy behavior in sync/update/rebuild scripts.
    """
    url = f"{_zotero_base()}/collections/{collection_key}/items/top?include=data"
    resp = requests.get(url, headers=_zotero_headers())
    resp.raise_for_status()

    keys: Set[str] = set()
    for item in resp.json():
        data = item.get("data", {}) or {}
        item_key = data.get("key")
        if item_key:
            keys.add(item_key)
    return keys


# --------- Notion helpers (stubs - Notion support removed) --------- #


def iter_literature_pages() -> List[Dict[str, Any]]:
    """Stub: Notion support removed. Returns empty list."""
    logger.debug("[MAINTENANCE][ZOTERO-UNIVERSE] iter_literature_pages() is a no-op (Notion removed)")
    return []


def delete_pinecone_vectors_for_item(item_key: str) -> None:
    """
    Stub: Pinecone support removed.
    
    Previously deleted all Pinecone vectors for a Zotero item.
    Use pgvector instead.
    """
    logger.debug(
        "[MAINTENANCE][ZOTERO-UNIVERSE] delete_pinecone_vectors_for_item() is a no-op (Pinecone deprecated)"
    )


def delete_rag_pages_for_lit_page(lit_page_id: str) -> None:
    """Stub: Notion support removed. No-op."""
    logger.debug("[MAINTENANCE][ZOTERO-UNIVERSE] delete_rag_pages_for_lit_page() is a no-op (Notion removed)")


def archive_literature_page(
    lit_page_id: str, title_text: str = "", item_key: str = ""
) -> None:
    """Stub: Notion support removed. No-op."""
    logger.debug("[MAINTENANCE][ZOTERO-UNIVERSE] archive_literature_page() is a no-op (Notion removed)")


# --------- cleanup_deleted_items.py --------- #


def cleanup_deleted_items() -> None:
    """
    Clean up deleted Zotero items:

    - Scan Literature DB.
    - For each page, read 'Zotero Item Key'.
    - If Zotero /items/{key} returns 404:
        - Delete Pinecone vectors.
        - Delete RAG Engine chunk pages.
        - Archive the Literature page.
    Mirrors cleanup_deleted_items.py.
    """
    pages = iter_literature_pages()
    kept = 0
    deleted = 0

    print("\n🔎 Cleaning up deleted Zotero items...\n")

    for page in pages:
        props = page.get("properties", {})
        z_prop = props.get("Zotero Item Key", {}).get("rich_text", [])
        if not z_prop:
            continue  # ignore non-Zotero pages

        item_key = z_prop[0].get("plain_text", "").strip()
        if not item_key:
            continue

        title_prop = props.get("Title", {}).get("title", [])
        title = title_prop[0]["plain_text"] if title_prop else "(untitled)"

        if zotero_item_exists(item_key):
            kept += 1
            continue

        print(
            f"❌ Zotero item missing for {title} ({item_key}) – deleting associated data..."
        )

        # 1) Pinecone
        delete_pinecone_vectors_for_item(item_key)

        # 2) RAG pages
        lit_page_id = page["id"]
        delete_rag_pages_for_lit_page(lit_page_id)

        # 3) Archive Literature page
        archive_literature_page(lit_page_id, title_text=title, item_key=item_key)

        deleted += 1

    print("\n=====================================")
    print("   ✅ Cleanup complete.")
    print(f"   Kept:    {kept} literature pages")
    print(f"   Deleted: {deleted} literature pages (and associated RAG+vectors)")
    print("=====================================\n")


# --------- sync_collection_state.py --------- #


def sync_collection_state(collection_key: str) -> None:
    """
    Enforce that ONLY items currently in the given Zotero collection exist in:
    - Literature DB
    - RAG Engine DB
    - Pinecone

    Mirrors sync_collection_state.py.
    """
    current_keys = collection_item_keys(collection_key)
    pages = iter_literature_pages()

    kept = 0
    deleted = 0

    print(
        f"\n🔄 Syncing collection-defined universe for collection {collection_key}...\n"
    )
    print(f"   Collection currently has {len(current_keys)} parent item(s).")

    for page in pages:
        props = page.get("properties", {})
        z_prop = props.get("Zotero Item Key", {}).get("rich_text", [])
        if not z_prop:
            continue

        item_key = z_prop[0].get("plain_text", "").strip()
        if not item_key:
            continue

        title_prop = props.get("Title", {}).get("title", [])
        title = title_prop[0]["plain_text"] if title_prop else "(untitled)"

        lit_page_id = page["id"]

        if item_key in current_keys and zotero_item_exists(item_key):
            kept += 1
            continue

        print(
            f"❌ Item {title} ({item_key}) not in collection or missing in Zotero → deleting."
        )

        delete_pinecone_vectors_for_item(item_key)
        delete_rag_pages_for_lit_page(lit_page_id)
        archive_literature_page(lit_page_id, title_text=title, item_key=item_key)

        deleted += 1

    print("\n=====================================")
    print("   ✅ Sync complete.")
    print(f"   Kept:    {kept} literature pages")
    print(f"   Deleted: {deleted} literature pages (and associated RAG+vectors)")
    print("=====================================\n")


# --------- update_collection_universe.py --------- #


def update_collection_universe(collection_key: str) -> None:
    """
    Incremental + delete-sync for a collection:

    1) Incrementally ingest all items in the collection
       (idempotent via ingest_zotero_item + md5 checks).
    2) Delete items that are no longer in the collection or missing in Zotero.

    Mirrors update_collection_universe.py.
    """
    print("\n📈 Incrementally ingesting collection items...")
    incremental_ingest_collection(
        collection_key=collection_key, parent_type="Literature"
    )

    print("\n🧹 Deleting items no longer in the collection or missing...")
    sync_collection_state(collection_key)


# --------- rebuild_collection_universe.py --------- #


def rebuild_collection_universe(collection_key: str) -> None:
    """
    Rebuild the entire RAG universe for a Zotero collection (Option A):

    1) For each item in the collection:
         - Delete its vectors and RAG pages.
         - Re-ingest via ingest_zotero_item() (via resync_collection).
    2) Delete items that are no longer in the collection or missing in Zotero.

    Mirrors rebuild_collection_universe.py.
    """

    print("\n🔁 Rebuilding collection universe (full resync)...")
    resync_collection(
        collection_key=collection_key,
        parent_type="Literature",
        force=True,
    )
    print("\n🧹 Cleaning up items no longer in collection or missing...")
    sync_collection_state(collection_key)
