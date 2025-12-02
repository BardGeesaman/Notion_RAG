# amprenta_rag/ingestion/zotero_api.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any, List, Optional

import requests

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class ZoteroItem:
    key: str
    title: str
    item_type: str
    abstract: str
    doi: Optional[str]
    url: Optional[str]
    date: Optional[str]
    journal: Optional[str]
    tags: List[str]


def _zotero_headers() -> Dict[str, str]:
    cfg = get_config().zotero
    return {
        "Zotero-API-Key": cfg.api_key,
        "Accept": "application/json",
    }


def _zotero_base() -> str:
    cfg = get_config().zotero
    return f"https://api.zotero.org/{cfg.library_type}s/{cfg.library_id}"


def fetch_zotero_item(item_key: str) -> ZoteroItem:
    """
    Fetch the main metadata for a Zotero item.
    """
    url = f"{_zotero_base()}/items/{item_key}?include=data"
    resp = requests.get(url, headers=_zotero_headers())
    resp.raise_for_status()
    data = resp.json().get("data", {})

    title = data.get("title") or "(untitled)"
    item_type = data.get("itemType") or "unknown"
    abstract = data.get("abstractNote") or ""
    doi = data.get("DOI") or None
    url_str = data.get("url") or None
    date = data.get("date") or None
    journal = data.get("publicationTitle") or data.get("journalAbbreviation") or None
    tags = [t.get("tag") for t in data.get("tags", []) if t.get("tag")]

    return ZoteroItem(
        key=item_key,
        title=title,
        item_type=item_type,
        abstract=abstract,
        doi=doi,
        url=url_str,
        date=date,
        journal=journal,
        tags=tags,
    )


def fetch_zotero_children(item_key: str) -> List[Dict[str, Any]]:
    """
    Fetch *all* children (attachments + notes) for a Zotero item.
    """
    url = f"{_zotero_base()}/items/{item_key}/children"
    resp = requests.get(url, headers=_zotero_headers())
    resp.raise_for_status()
    return [c.get("data", {}) for c in resp.json()]


def fetch_zotero_attachments(item_key: str) -> List[Dict[str, Any]]:
    """
    All child items of type 'attachment'.
    """
    return [d for d in fetch_zotero_children(item_key) if d.get("itemType") == "attachment"]


def fetch_zotero_notes(item_key: str) -> List[Dict[str, Any]]:
    """
    All child items of type 'note'.
    """
    return [d for d in fetch_zotero_children(item_key) if d.get("itemType") == "note"]


def download_zotero_file(att_key: str) -> bytes:
    """
    Download raw bytes for an attachment from Zotero.
    """
    url = f"{_zotero_base()}/items/{att_key}/file"
    cfg = get_config().zotero
    resp = requests.get(url, headers={"Zotero-API-Key": cfg.api_key})
    resp.raise_for_status()
    return resp.content