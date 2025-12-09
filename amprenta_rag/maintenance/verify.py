# amprenta_rag/maintenance/verify.py
"""
Verification utilities for maintenance operations.

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from typing import Any, Dict, List

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _query_db(db_id: str, payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Stub: Notion support removed. Returns empty list."""
    logger.debug("[MAINTENANCE][VERIFY] _query_db() is a no-op (Notion removed)")
    return []


def _count_rag_chunks_for_lit_page(lit_page_id: str) -> int:
    cfg = get_config().notion
    rag_db_id = cfg.rag_db_id
    if not rag_db_id:
        return 0

    payload = {
        "page_size": 100,
        "filter": {
            "property": "Parent Item",
            "relation": {"contains": lit_page_id},
        },
    }
    chunks = _query_db(rag_db_id, payload)
    return len(chunks)


def _count_rag_chunks_for_email_page(email_page_id: str) -> int:
    cfg = get_config().notion
    rag_db_id = cfg.rag_db_id
    if not rag_db_id:
        return 0

    payload = {
        "page_size": 100,
        "filter": {
            "property": "Parent Email/Note Item",
            "relation": {"contains": email_page_id},
        },
    }
    chunks = _query_db(rag_db_id, payload)
    return len(chunks)


def verify_rag_metadata() -> None:
    """
    Verify RAG metadata consistency between Notion Literature/Email DBs and RAG DB.

    Mirrors rag_verify_metadata.py.
    """
    cfg = get_config().notion
    lit_db_id = cfg.lit_db_id
    email_db_id = cfg.email_db_id

    print("\n🔎 Verifying RAG metadata consistency...\n")

    # 1) Literature pages
    if lit_db_id:
        lit_pages = _query_db(lit_db_id, {"page_size": 100})
        print("📚 Literature pages:")
        for page in lit_pages:
            props = page.get("properties", {})
            z_prop = props.get("Zotero Item Key", {}).get("rich_text", [])
            if not z_prop:
                continue
            item_key = z_prop[0].get("plain_text", "").strip()
            title_prop = props.get("Title", {}).get("title", [])
            title = title_prop[0]["plain_text"] if title_prop else "(untitled)"

            chunk_count = _count_rag_chunks_for_lit_page(page["id"])
            print(f"• [{item_key}] {title} -> {chunk_count} RAG chunks")

    # 2) Email/Note pages
    if email_db_id:
        payload = {
            "page_size": 100,
            "filter": {
                "property": "Embedding Status",
                "select": {"equals": "Embedded"},
            },
        }
        email_pages = _query_db(email_db_id, payload)
        print("\n📧 Email/Note pages (Embedded):")
        for page in email_pages:
            props = page.get("properties", {})
            title_prop = props.get("Title", {}).get("title", [])
            title = title_prop[0]["plain_text"] if title_prop else "(untitled)"
            type_prop = props.get("Type", {}).get("select", {})
            parent_type = type_prop.get("name", "Email")

            chunk_count = _count_rag_chunks_for_email_page(page["id"])
            print(f"• [{parent_type}] {title} -> {chunk_count} RAG chunks")

    print("\n✅ Verification complete.\n")
