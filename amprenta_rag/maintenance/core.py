# amprenta_rag/maintenance/core.py

from __future__ import annotations

import json
import time
from typing import Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[MAINTENANCE][CORE] notion_headers() deprecated - Notion support removed")
    return {}


def delete_all_pinecone_vectors() -> None:
    """
    Delete ALL vectors from the configured Pinecone index.
    Mirrors legacy delete_pinecone_data.py.
    """
    cfg = get_config()
    index = get_pinecone_index()
    index.delete(delete_all=True, namespace=cfg.pinecone.namespace)
    print("🧹 Pinecone vectors deleted (namespace:", cfg.pinecone.namespace, ")")


def archive_all_pages_in_db(db_id: str, label: str = "") -> None:
    """
    Archive all pages in a Notion database by setting archived=true.
    Mirrors the behavior used in clear_rag_engine.py and reset_all.py.
    """
    base_url = get_config().notion.base_url
    url = f"{base_url}/databases/{db_id}/query"

    payload = {"page_size": 100}
    total_archived = 0
    next_cursor: Optional[str] = None

    print(f"\n🧹 Archiving all pages in {label or 'Notion DB'} ({db_id})...")

    while True:
        body = payload.copy()
        if next_cursor:
            body["start_cursor"] = next_cursor

        resp = requests.post(url, headers=notion_headers(), data=json.dumps(body))
        resp.raise_for_status()
        data = resp.json()

        results = data.get("results", [])
        if not results:
            break

        for page in results:
            page_id = page["id"]
            upd = {"archived": True}
            upd_url = f"{base_url}/pages/{page_id}"

            resp2 = requests.patch(
                upd_url, headers=notion_headers(), data=json.dumps(upd)
            )
            if resp2.status_code >= 300:
                print("⚠️ Error archiving page:", resp2.text)
            else:
                total_archived += 1

            time.sleep(0.10)  # gentle rate limit

        if not data.get("has_more"):
            break
        next_cursor = data.get("next_cursor")

    print(f"✅ Done. Archived {total_archived} pages from {label or 'DB'}.")


def clear_rag_engine_db() -> None:
    """
    Archive all pages in the RAG Engine DB.
    Mirrors clear_rag_engine.py.
    """
    cfg = get_config().notion
    rag_db_id = cfg.rag_db_id
    if not rag_db_id:
        raise RuntimeError("NOTION_RAG_DB_ID is not configured in config.py")
    archive_all_pages_in_db(rag_db_id, label="RAG Engine")


def reset_all() -> None:
    """
    Fully reset the RAG system:

    1. Delete all vectors from Pinecone.
    2. Archive all pages in the RAG Engine DB.
    3. Archive all pages in the Literature DB.

    Mirrors reset_all.py.
    """
    cfg = get_config()
    print("\n🚨 FULL RAG RESET STARTING...\n")

    # 1) Pinecone
    delete_all_pinecone_vectors()

    # 2) RAG Engine DB
    if cfg.notion.rag_db_id:
        archive_all_pages_in_db(cfg.notion.rag_db_id, label="RAG Engine")
    else:
        print("⚠️ NOTION_RAG_DB_ID not set; skipping RAG Engine DB archive.")

    # 3) Literature DB
    if cfg.notion.lit_db_id:
        archive_all_pages_in_db(cfg.notion.lit_db_id, label="Literature")
    else:
        print("⚠️ NOTION_LIT_DB_ID not set; skipping Literature DB archive.")

    print("\n🎉 RAG reset complete. You now have a clean slate for re-ingestion.\n")
