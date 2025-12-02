# amprenta_rag/ingestion/__init__.py

from .zotero_ingest import ingest_zotero_item
from .zotero_collection import incremental_ingest_collection, resync_collection
from .email_ingestion import (
    ingest_email,
    batch_ingest_emails,
    delete_email_and_chunks,
    cleanup_orphaned_chunks,
)

__all__ = [
    "ingest_zotero_item",
    "incremental_ingest_collection",
    "resync_collection",
    "ingest_email",
    "batch_ingest_emails",
    "delete_email_and_chunks",
    "cleanup_orphaned_chunks",
]