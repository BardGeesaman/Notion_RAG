"""
Zotero item ingestion pipeline.

This package provides functions for ingesting Zotero items into Notion and Pinecone,
including attachment processing, note processing, feature extraction, and signature detection.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.zotero.ingestion import ingest_zotero_item

__all__ = [
    "ingest_zotero_item",
]

