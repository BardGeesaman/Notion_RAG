"""
Compatibility wrapper for the refactored Zotero ingestion pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from zotero_ingest.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.zotero import ingest_zotero_item

__all__ = [
    "ingest_zotero_item",
]
