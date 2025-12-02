# amprenta_rag/ingestion/notion_pages.py

"""
Notion API helpers for page creation and management.

This module provides functions for:
- Creating and updating Literature DB pages
- Creating RAG chunk pages
- Fetching and extracting content from Notion pages
- Email/Note specific Notion operations
"""

from __future__ import annotations

from typing import Dict, Any, Optional, List

import json
import re
import time
from datetime import datetime, timezone

import requests

from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.ingestion.zotero_api import ZoteroItem
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Rate limiting for Notion API
REQUESTS_PER_SECOND = 2.0


def _require_notion_db_id(name: str, value: Optional[str]) -> str:
    if not value:
        raise RuntimeError(f"Notion {name} DB ID is not configured in config.py")
    return value


def ensure_literature_page(item: ZoteroItem, parent_type: str) -> str:
    """
    Find or create a Literature DB page for this Zotero item.

    Matches your schema:
      - Title (title)
      - Source Type (select)  -> "Paper"
      - Embedding Status (select) -> "Not Embedded"
      - Zotero Item Key (rich_text)
      - Abstract, Journal, DOI (url), Year, URL, Tags (optional)
    """
    cfg = get_config().notion
    lit_db_id = _require_notion_db_id("LIT", cfg.lit_db_id)

    # 1) Try find by Zotero Item Key
    query_url = f"{cfg.base_url}/databases/{lit_db_id}/query"
    payload = {
        "page_size": 1,
        "filter": {
            "property": "Zotero Item Key",
            "rich_text": {"equals": item.key},
        },
    }

    try:
        resp = requests.post(query_url, headers=notion_headers(), json=payload)
        if resp.status_code >= 300:
            logger.warning(
                "[NOTION] Literature query error for item %s: %s %s",
                item.key,
                resp.status_code,
                resp.text,
            )
            resp.raise_for_status()
    except Exception as e:
        logger.error(
            "[NOTION] Notion API error querying literature page for item %s: %r",
            item.key,
            e,
        )
        raise

    results = resp.json().get("results", [])
    if results:
        return results[0]["id"]

    # 2) Create new
    logger.info("[NOTION] Creating new Literature page for Zotero item %s", item.key)

    year_val = None
    if item.date:
        m = re.search(r"(19|20)\d{2}", item.date)
        if m:
            year_val = int(m.group(0))

    props: Dict[str, Any] = {
        "Title": {"title": [{"text": {"content": item.title or "Untitled"}}]},
        "Source Type": {"select": {"name": "Paper"}},
        "Embedding Status": {"select": {"name": "Not Embedded"}},
        "Zotero Item Key": {"rich_text": [{"text": {"content": item.key}}]},
    }

    if item.abstract:
        props["Abstract"] = {
            "rich_text": [{"text": {"content": item.abstract}}],
        }
    if item.journal:
        props["Journal"] = {
            "rich_text": [{"text": {"content": item.journal}}],
        }
    if item.doi:
        doi_str = item.doi.strip()
        if doi_str and not doi_str.startswith("http"):
            doi_str = "https://doi.org/" + doi_str
        props["DOI"] = {"url": doi_str}
    if year_val is not None:
        props["Year"] = {"number": year_val}
    if item.url:
        props["URL"] = {"url": item.url}
    if item.tags:
        props["Tags"] = {"multi_select": [{"name": t} for t in item.tags]}

    create_payload = {
        "parent": {"database_id": lit_db_id},
        "properties": props,
    }

    try:
        resp = requests.post(
            f"{cfg.base_url}/pages",
            headers=notion_headers(),
            json=create_payload,
        )
        if resp.status_code >= 300:
            logger.error(
                "[NOTION] Literature create error for item %s: %s %s",
                item.key,
                resp.status_code,
                resp.text,
            )
            resp.raise_for_status()
    except Exception as e:
        logger.error(
            "[NOTION] Notion API error creating literature page for item %s: %r",
            item.key,
            e,
        )
        raise

    return resp.json()["id"]


def update_literature_page(parent_page_id: str, when_iso: Optional[str] = None) -> None:
    """
    Set Embedding Status = 'Embedded' and update Last Embedded.
    """
    cfg = get_config().notion
    if when_iso is None:
        when_iso = datetime.now(timezone.utc).isoformat()

    payload = {
        "properties": {
            "Embedding Status": {"select": {"name": "Embedded"}},
            "Last Embedded": {"date": {"start": when_iso}},
        },
    }

    try:
        resp = requests.patch(
            f"{cfg.base_url}/pages/{parent_page_id}",
            headers=notion_headers(),
            data=json.dumps(payload),
        )
        if resp.status_code >= 300:
            logger.warning("[NOTION] Literature update error for page %s: %s", parent_page_id, resp.text)
    except Exception as e:
        logger.error(
            "[NOTION] Notion API error updating literature page %s: %r",
            parent_page_id,
            e,
        )
        raise


def create_rag_chunk_page(
    chunk_id: str,
    chunk_text: str,
    parent_type: str,
    parent_page_id: str,
    order: int,
    when_iso: str,
    parent_relation_key: Optional[str] = None,
) -> str:
    """
    Create a chunk page in the RAG DB, matching your original schema:

      - Chunk ID (title)
      - Chunk Text (rich_text)
      - Parent Type (select)
      - Parent Item (relation) or Parent Email/Note Item (relation) for Email/Note
      - Order (number)
      - Embedding Vector ID (rich_text)
      - Last Embedded (date)

    Args:
        chunk_id: Unique identifier for the chunk
        chunk_text: Text content of the chunk (truncated to 1900 chars)
        parent_type: Type of parent (e.g., "Literature", "Email", "Note")
        parent_page_id: Notion page ID of the parent (without dashes)
        order: Order/index of this chunk within the parent
        when_iso: ISO timestamp for "Last Embedded"
        parent_relation_key: Optional relation property name. If None, defaults to
            "Parent Email/Note Item" for Email/Note types, "Parent Item" otherwise.
    """
    cfg = get_config().notion
    rag_db_id = _require_notion_db_id("RAG", cfg.rag_db_id)

    # Determine relation property based on parent_type if not explicitly provided
    if parent_relation_key is None:
        if parent_type in ("Email", "Note"):
            parent_relation_key = "Parent Email/Note Item"
        else:
            parent_relation_key = "Parent Item"

    props: Dict[str, Any] = {
        "Chunk ID": {"title": [{"text": {"content": chunk_id}}]},
        "Chunk Text": {"rich_text": [{"text": {"content": chunk_text[:1900]}}]},
        "Parent Type": {"select": {"name": parent_type}},
        parent_relation_key: {"relation": [{"id": parent_page_id}]},
        "Order": {"number": order},
        "Embedding Vector ID": {"rich_text": [{"text": {"content": chunk_id}}],
        },
        "Last Embedded": {"date": {"start": when_iso}},
    }

    payload = {
        "parent": {"database_id": rag_db_id},
        "properties": props,
    }

    try:
        resp = requests.post(
            f"{cfg.base_url}/pages",
            headers=notion_headers(),
            data=json.dumps(payload),
        )

        if resp.status_code >= 300:
            logger.warning(
                "[NOTION] RAG chunk create error for chunk %s: %s %s",
                chunk_id,
                resp.status_code,
                resp.text,
            )
            resp.raise_for_status()
    except Exception as e:
        logger.error(
            "[NOTION] Notion API error creating RAG chunk page %s: %r",
            chunk_id,
            e,
        )
        raise

    page_id = resp.json()["id"].replace("-", "")
    return page_id


# ----------------- Email/Note specific Notion helpers ----------------- #


def fetch_not_embedded_emails() -> List[Dict[str, Any]]:
    """
    Fetch all email/note pages from the Email DB whose 'Embedding Status' is 'Not Embedded'.
    
    Returns:
        List of Notion page objects from the Email DB
    """
    cfg = get_config().notion
    if not cfg.email_db_id:
        raise RuntimeError("Notion Email DB ID is not configured in config.py")
    
    db_id = cfg.email_db_id
    query_url = f"{cfg.base_url}/databases/{db_id}/query"

    logger.info("[NOTION] Fetching emails with 'Not Embedded' status...")

    payload: Dict[str, Any] = {
        "page_size": 100,
        "filter": {
            "property": "Embedding Status",
            "select": {"equals": "Not Embedded"},
        },
    }

    results: List[Dict[str, Any]] = []
    next_cursor: Optional[str] = None

    while True:
        body = payload.copy()
        if next_cursor:
            body["start_cursor"] = next_cursor

        try:
            resp = requests.post(query_url, headers=notion_headers(), json=body)
            if resp.status_code >= 300:
                logger.error("[NOTION] Error querying Email DB: %s", resp.text)
                resp.raise_for_status()
        except Exception as e:
            logger.error("[NOTION] Notion API error querying Email DB: %r", e)
            raise

        data = resp.json()
        batch = data.get("results", [])
        results.extend(batch)

        if not data.get("has_more"):
            break
        next_cursor = data.get("next_cursor")

        time.sleep(1.0 / REQUESTS_PER_SECOND)

    logger.info("[NOTION] Got %d not-embedded emails", len(results))
    return results


def extract_page_content(page_id: str) -> str:
    """
    Fetch the full page content (blocks) from Notion and extract text.
    
    Extracts text from various Notion block types (paragraphs, headings, lists, etc.)
    and returns a plain text representation. This is email/note-specific as it
    operates on Notion page blocks rather than file attachments.
    
    Args:
        page_id: Notion page ID (without dashes) to extract content from
        
    Returns:
        Plain text representation of the page content
    """
    cfg = get_config().notion
    url = f"{cfg.base_url}/blocks/{page_id}/children"
    texts: List[str] = []

    next_cursor: Optional[str] = None
    while True:
        params = {}
        if next_cursor:
            params["start_cursor"] = next_cursor

        try:
            resp = requests.get(url, headers=notion_headers(), params=params or None)
            resp.raise_for_status()
        except Exception as e:
            logger.error(
                "[NOTION] Notion API error fetching page blocks for page %s: %r",
                page_id,
                e,
            )
            raise
        data = resp.json()

        blocks = data.get("results", [])
        for block in blocks:
            block_type = block.get("type")

            def _plain(key: str) -> str:
                rich = block.get(key, {}).get("rich_text", []) or []
                return "".join(t.get("plain_text", "") for t in rich)

            if block_type == "paragraph":
                texts.append(_plain("paragraph"))
            elif block_type == "heading_1":
                texts.append("# " + _plain("heading_1"))
            elif block_type == "heading_2":
                texts.append("## " + _plain("heading_2"))
            elif block_type == "heading_3":
                texts.append("### " + _plain("heading_3"))
            elif block_type == "bulleted_list_item":
                texts.append("• " + _plain("bulleted_list_item"))
            elif block_type == "numbered_list_item":
                texts.append(_plain("numbered_list_item"))
            elif block_type == "quote":
                texts.append(_plain("quote"))
            elif block_type == "code":
                texts.append(_plain("code"))

        if not data.get("has_more"):
            break
        next_cursor = data.get("next_cursor")

        time.sleep(1.0 / REQUESTS_PER_SECOND)

    return "\n".join(texts)


def update_email_page(parent_page_id: str, when_iso: str) -> None:
    """
    Update email/note page status after successful ingestion.
    
    Sets Embedding Status to 'Embedded' and updates Last Embedded timestamp.
    This marks the email/note as successfully processed in the RAG system.
    
    Args:
        parent_page_id: Notion page ID of the email/note (without dashes)
        when_iso: ISO timestamp string for when embedding occurred
    """
    cfg = get_config().notion
    payload = {
        "properties": {
            "Embedding Status": {"select": {"name": "Embedded"}},
            "Last Embedded": {"date": {"start": when_iso}},
        },
    }

    try:
        resp = requests.patch(
            f"{cfg.base_url}/pages/{parent_page_id}",
            headers=notion_headers(),
            data=json.dumps(payload),
        )

        if resp.status_code >= 300:
            logger.error(
                "[NOTION] Email update error for page %s: %s",
                parent_page_id,
                resp.text,
            )
            raise Exception(f"Failed to update email status: {resp.text}")
    except Exception as e:
        logger.error(
            "[NOTION] Notion API error updating email page %s: %r",
            parent_page_id,
            e,
        )
        raise