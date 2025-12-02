# amprenta_rag/ingestion/pinecone_utils.py

"""
Pinecone utility functions.

This module provides:
- Metadata sanitization for Pinecone compatibility
- Idempotency checks (attachment/note already ingested)
- Query helpers for Pinecone operations
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional

from amprenta_rag.config import get_config
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def sanitize_metadata(meta: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ensure all metadata values are allowed by Pinecone:
      - string
      - number (int/float)
      - boolean
      - list of strings

    Drop keys with None or empty lists.
    Convert non-string, non-primitive list items to strings.
    """
    cleaned: Dict[str, Any] = {}

    for key, value in meta.items():
        if value is None:
            # Pinecone does not allow null metadata
            continue

        # Primitive types are fine
        if isinstance(value, (str, bool, int, float)):
            cleaned[key] = value
            continue

        # Lists need to be list of strings
        if isinstance(value, list):
            string_list = [str(v) for v in value if v is not None]
            if string_list:
                cleaned[key] = string_list
            continue

        # Anything else (e.g. dicts) – convert to string for now
        cleaned[key] = str(value)

    return cleaned


def _query_pinecone_by_filter(filter_obj: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Helper: query Pinecone by metadata filter using a dummy vector of
    the right dimension. We rely on the index being 3072-dim (current index).
    """
    index = get_pinecone_index()
    cfg = get_config()

    dummy = [0.0] * 3072
    try:
        res = index.query(
            vector=dummy,
            top_k=10,
            include_metadata=True,
            namespace=cfg.pinecone.namespace,
            filter=filter_obj,
        )
    except Exception as e:
        logger.error("[PINECONE] Pinecone API error querying by filter: %r", e)
        raise

    matches = getattr(res, "matches", None)
    if matches is None:
        matches = res.get("matches", [])
    return list(matches)


def attachment_already_ingested(
    item_key: str,
    attachment_key: str,
    current_md5: Optional[str] = None,
) -> bool:
    """
    Return True if this (item_key, attachment_key) is already ingested AND
    md5 matches. If md5 changed, delete old vectors and return False.
    """
    matches = _query_pinecone_by_filter(
        {
            "zotero_item_key": {"$eq": item_key},
            "attachment_key": {"$eq": attachment_key},
        }
    )
    if not matches:
        return False

    if not current_md5:
        return True

    m0 = matches[0]
    meta = getattr(m0, "metadata", None) or m0.get("metadata", {})
    stored_md5 = meta.get("attachment_md5")

    index = get_pinecone_index()
    cfg = get_config()

    if (not stored_md5) or (stored_md5 != current_md5):
        logger.warning(
            "[PINECONE] Attachment %s for item %s has changed (stored md5=%r, current md5=%r); "
            "deleting old vectors and re-ingesting.",
            attachment_key,
            item_key,
            stored_md5,
            current_md5,
        )
        try:
            index.delete(
                filter={
                    "zotero_item_key": {"$eq": item_key},
                    "attachment_key": {"$eq": attachment_key},
                },
                namespace=cfg.pinecone.namespace,
            )
        except Exception as e:
            logger.error(
                "[PINECONE] Pinecone API error deleting old vectors for attachment %s (item %s): %r",
                attachment_key,
                item_key,
                e,
            )
            raise
        return False

    return True


def note_already_ingested(
    item_key: str,
    note_key: str,
    note_hash: str,
) -> bool:
    """
    Same pattern as attachments: check if note with this hash is already in Pinecone.
    If hash changed, delete old vectors.
    """
    matches = _query_pinecone_by_filter(
        {
            "zotero_item_key": {"$eq": item_key},
            "note_key": {"$eq": note_key},
        }
    )
    if not matches:
        return False

    m0 = matches[0]
    meta = getattr(m0, "metadata", None) or m0.get("metadata", {})
    stored_hash = meta.get("note_hash")

    index = get_pinecone_index()
    cfg = get_config()

    if (not stored_hash) or (stored_hash != note_hash):
        logger.warning(
            "[PINECONE] Note %s for item %s has changed (stored=%r, current=%r); deleting old vectors and re-ingesting.",
            note_key,
            item_key,
            stored_hash,
            note_hash,
        )
        try:
            index.delete(
                filter={
                    "zotero_item_key": {"$eq": item_key},
                    "note_key": {"$eq": note_key},
                },
                namespace=cfg.pinecone.namespace,
            )
        except Exception as e:
            logger.error(
                "[PINECONE] Pinecone API error deleting old vectors for note %s (item %s): %r",
                note_key,
                item_key,
                e,
            )
            raise
        return False

    return True