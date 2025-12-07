"""
Postgres integration for RAG chunk creation.

Provides functions to create RAG chunks in Postgres,
replacing Notion-based chunk page creation.
"""

from __future__ import annotations

from typing import Any, Dict, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import RAGChunk
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_rag_chunk_in_postgres(
    chunk_id: str,
    chunk_text: str,
    source_type: str,
    source_id: Optional[UUID] = None,
    source_name: Optional[str] = None,
    chunk_index: int = 0,
    snippet: Optional[str] = None,
    zotero_item_key: Optional[str] = None,
    note_key: Optional[str] = None,
    attachment_key: Optional[str] = None,
    note_hash: Optional[str] = None,
    attachment_hash: Optional[str] = None,
    chunk_metadata: Optional[Dict[str, Any]] = None,
    literature_id: Optional[UUID] = None,
    email_id: Optional[UUID] = None,
    notion_page_id: Optional[str] = None,
    db: Optional[Session] = None,
) -> UUID:
    """
    Create a RAG chunk record in Postgres.

    Args:
        chunk_id: Unique chunk identifier (original format)
        chunk_text: Full text content of the chunk
        source_type: Type of source (Literature, Dataset, Email, etc.)
        source_id: UUID of the source record (literature, dataset, email, etc.)
        source_name: Name of the source
        chunk_index: Index of chunk within source (0-based)
        snippet: Short preview text
        zotero_item_key: Zotero item key (if from Zotero)
        note_key: Zotero note key (if from note)
        attachment_key: Zotero attachment key (if from attachment)
        note_hash: MD5 hash of note content
        attachment_hash: MD5 hash of attachment content
        chunk_metadata: Additional metadata as JSON dict
        literature_id: UUID of Literature record (if applicable)
        email_id: UUID of Email record (if applicable)
        notion_page_id: Optional Notion page ID for backward compatibility
        db: Optional database session

    Returns:
        UUID of the created RAG chunk record

    Example:
        >>> chunk_id = create_rag_chunk_in_postgres(
        ...     chunk_id="abc123_note_001",
        ...     chunk_text="This is a chunk of text...",
        ...     source_type="Literature",
        ...     source_name="Paper Title",
        ...     chunk_index=0,
        ... )
        >>> isinstance(chunk_id, UUID)
        True
    """
    if db is None:
        with get_db() as session:
            return create_rag_chunk_in_postgres(
                chunk_id=chunk_id,
                chunk_text=chunk_text,
                source_type=source_type,
                source_id=source_id,
                source_name=source_name,
                chunk_index=chunk_index,
                snippet=snippet,
                zotero_item_key=zotero_item_key,
                note_key=note_key,
                attachment_key=attachment_key,
                note_hash=note_hash,
                attachment_hash=attachment_hash,
                chunk_metadata=chunk_metadata,
                literature_id=literature_id,
                email_id=email_id,
                notion_page_id=notion_page_id,
                db=session,
            )

    # Check if chunk already exists
    existing = db.query(RAGChunk).filter(RAGChunk.chunk_id == chunk_id).first()
    if existing:
        logger.debug(
            "[INGEST][RAG-CHUNK][POSTGRES] Chunk %s already exists: %s",
            chunk_id,
            existing.id,
        )
        return existing.id

    # Create new chunk
    chunk = RAGChunk(
        chunk_id=chunk_id,
        chunk_text=chunk_text,
        chunk_index=chunk_index,
        snippet=snippet,
        source_type=source_type,
        source_id=source_id,
        source_name=source_name,
        zotero_item_key=zotero_item_key,
        note_key=note_key,
        attachment_key=attachment_key,
        note_hash=note_hash,
        attachment_hash=attachment_hash,
        chunk_metadata=chunk_metadata,
        literature_id=literature_id,
        email_id=email_id,
        notion_page_id=notion_page_id,
    )

    db.add(chunk)
    db.commit()
    db.refresh(chunk)

    logger.debug(
        "[INGEST][RAG-CHUNK][POSTGRES] Created chunk %s (ID: %s)",
        chunk_id,
        chunk.id,
    )

    return chunk.id
