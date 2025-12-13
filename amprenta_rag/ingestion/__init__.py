# amprenta_rag/ingestion/__init__.py

from .email_ingestion import (
    batch_ingest_emails,
    cleanup_orphaned_chunks,
    delete_email_and_chunks,
    ingest_email,
)
from .experiments_ingestion import ingest_experiment
from .zotero_collection import incremental_ingest_collection, resync_collection
from .zotero_ingest import ingest_zotero_item

def ingest_dataset(*args, **kwargs):
    """
    Backwards-compatible dataset ingestion entrypoint.

    The legacy `dataset_ingestion.py` module was removed as part of the Postgres-only
    migration. This function lazily delegates to the Postgres ingestion pipeline.
    """
    from .postgres_dataset_ingestion import ingest_dataset_from_postgres

    return ingest_dataset_from_postgres(*args, **kwargs)

__all__ = [
    "ingest_zotero_item",
    "incremental_ingest_collection",
    "resync_collection",
    "ingest_email",
    "batch_ingest_emails",
    "delete_email_and_chunks",
    "cleanup_orphaned_chunks",
    "ingest_experiment",
    "ingest_dataset",
]
