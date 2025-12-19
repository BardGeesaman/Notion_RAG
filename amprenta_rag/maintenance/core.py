# amprenta_rag/maintenance/core.py
"""
Core maintenance utilities.

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations


from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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
    Stub: Notion support removed. No-op.

    Previously archived all pages in a Notion database.
    """
    logger.debug("[MAINTENANCE][CORE] archive_all_pages_in_db() is a no-op (Notion removed)")
    print("⚠️ archive_all_pages_in_db() is a no-op (Notion support removed)")


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
