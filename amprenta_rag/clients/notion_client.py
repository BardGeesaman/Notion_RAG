# amprenta_rag/clients/notion_client.py

from __future__ import annotations

from typing import Any, Dict, List

import requests

from amprenta_rag.config import get_config


def notion_headers() -> Dict[str, str]:
    cfg = get_config().notion
    return {
        "Authorization": f"Bearer {cfg.api_key}",
        "Notion-Version": cfg.version,
        "Content-Type": "application/json",
    }


def _fetch_chunk_text_property(page_id: str) -> str:
    """
    Old behavior: read the 'Chunk Text' rich_text property from the page.
    This matches your original rag_query.py.
    """
    cfg = get_config().notion
    url = f"{cfg.base_url}/pages/{page_id}"
    resp = requests.get(url, headers=notion_headers())
    if resp.status_code >= 300:
        # Mirror your old behavior: just return "" on error
        return ""
    page = resp.json()
    props = page.get("properties", {})
    chunk_prop = props.get("Chunk Text", {})
    rich = chunk_prop.get("rich_text", []) or []
    return "".join(r.get("plain_text", "") for r in rich)


def _fetch_block_text(page_id: str) -> str:
    """
    Fallback: read all child blocks and concatenate their plain text.
    Useful if future chunks are stored as page content instead of a property.
    """
    cfg = get_config().notion
    url = f"{cfg.base_url}/blocks/{page_id}/children"
    texts: List[str] = []

    while url:
        resp = requests.get(url, headers=notion_headers())
        if resp.status_code >= 300:
            break
        data = resp.json()
        for block in data.get("results", []):
            bt = _block_plain_text(block)
            if bt:
                texts.append(bt)
        if data.get("has_more") and data.get("next_cursor"):
            url = f"{cfg.base_url}/blocks/{page_id}/children?start_cursor={data['next_cursor']}"
        else:
            url = None

    return "\n".join(texts).strip()


def _block_plain_text(block: Dict[str, Any]) -> str:
    """
    Extract plain text from a generic Notion block.
    Currently only handles paragraph blocks, which is fine for RAG chunks.
    """
    block_type = block.get("type")
    if not block_type:
        return ""
    rich = block.get(block_type, {}).get("rich_text") or []
    return "".join(span.get("plain_text", "") for span in rich)


def get_page_text(page_id: str) -> str:
    """
    Unified helper used by the RAG query engine:

    1. First try the 'Chunk Text' property (your current schema).
    2. If empty, fall back to concatenating child block text.
    """
    # 1) Old behavior – property-based chunks
    text = _fetch_chunk_text_property(page_id)
    if text.strip():
        return text.strip()

    # 2) Fallback – blocks-based content
    return _fetch_block_text(page_id)
