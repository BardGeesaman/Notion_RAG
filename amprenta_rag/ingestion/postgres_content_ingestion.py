"""
Postgres-only content ingestion module.

Handles ingestion of various content types (emails, Zotero/literature)
directly to Pinecone without Notion dependencies.

For content that doesn't need to be stored in Postgres (emails, literature),
this module provides direct-to-Pinecone ingestion for fast, scalable operation.
"""

from __future__ import annotations

import hashlib
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import uuid4

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import extract_features_from_text
from amprenta_rag.utils.metadata import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _calculate_content_hash(content: str) -> str:
    """Calculate SHA256 hash of content for change detection."""
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _content_already_ingested(
    content_id: str,
    force: bool = False,
) -> tuple[bool, Optional[str]]:
    """
    Check if content with this content_id is already ingested.

    Returns:
        Tuple of (is_already_ingested, stored_hash)
        If content changed, returns (False, stored_hash) and caller should delete old vectors
    """
    if force:
        return False, None

    if not content_id:
        return False, None

    # NOTE: Removed pinecone-specific deduplication check
    # Content will be reprocessed if ingested multiple times
    # NOTE: pgvector-based deduplication tracked in ROADMAP
    return False, None


def _delete_content_vectors(content_id: str) -> None:
    """Delete all vectors for a given content_id from Pinecone."""
    cfg = get_config()
    store = get_vector_store()

    try:
        store.delete(
            filter={
                "content_id": {"$eq": content_id},
            },
            namespace=cfg.pinecone.namespace,
        )
        logger.info(
            "[INGEST][POSTGRES-CONTENT] Deleted old vectors for content %s",
            content_id[:8],
        )
    except Exception as e:
        logger.warning(
            "[INGEST][POSTGRES-CONTENT] Error deleting old vectors for content %s: %r",
            content_id[:8],
            e,
        )
        # Don't raise - continue with ingestion


def ingest_content_direct_to_pinecone(
    content: str,
    source_type: str,
    title: str,
    metadata: Optional[Dict[str, Any]] = None,
    content_id: Optional[str] = None,
    attachment_paths: Optional[List[Path]] = None,
    force: bool = False,
) -> List[str]:
    """
    Ingest content directly to Pinecone without Notion or Postgres storage.

    This function is designed for transient content (emails, literature) that
    doesn't need to be stored in Postgres. Content goes directly to Pinecone
    for fast, scalable ingestion.

    Idempotent: Only ingests if content is new or changed. Uses content_hash
    to detect changes. Existing content with same hash is skipped.

    Args:
        content: Full text content to ingest
        source_type: Type of content ("email", "literature", "note", etc.)
        title: Title of the content
        metadata: Optional metadata dictionary
        content_id: Optional unique identifier (if None, generates UUID)
        attachment_paths: Optional list of attachment file paths
        force: If True, re-ingest even if already ingested

    Returns:
        List of embedding/chunk IDs created (empty if skipped)

    Raises:
        Exception: If ingestion fails
    """
    if not content or len(content.strip()) < 50:
        logger.warning(
            "[INGEST][POSTGRES-CONTENT] Content too short (%d chars), skipping",
            len(content) if content else 0,
        )
        return []

    # Generate content ID if not provided
    if not content_id:
        content_id = str(uuid4()).replace("-", "")

    # Calculate content hash for change detection
    content_hash = _calculate_content_hash(content)

    # Check if already ingested (idempotency check)
    already_ingested, stored_hash = _content_already_ingested(
        content_id=content_id,
        current_hash=content_hash,
        force=force,
    )

    if already_ingested:
        logger.info(
            "[INGEST][POSTGRES-CONTENT] Content %s already ingested and unchanged; skipping",
            content_id[:8],
        )
        return []

    # Content is new or changed - delete old vectors if changed
    if stored_hash:
        _delete_content_vectors(content_id)

    logger.info(
        "[INGEST][POSTGRES-CONTENT] Ingesting %s content: %s (ID: %s)",
        source_type,
        title[:60],
        content_id[:8],
    )

    # Build base metadata
    base_meta: Dict[str, Any] = {
        "source": source_type.title(),
        "source_type": source_type,
        "content_id": content_id,
        "content_hash": content_hash,  # Store hash for change detection
        "title": title[:200],
        **(metadata or {}),
    }

    # Chunk and embed
    chunks = chunk_text(content)
    if not chunks:
        logger.warning(
            "[INGEST][POSTGRES-CONTENT] No chunks produced for content %s",
            content_id[:8],
        )
        return []

    logger.info(
        "[INGEST][POSTGRES-CONTENT] Generated %d chunk(s) for content %s",
        len(chunks),
        content_id[:8],
    )

    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES-CONTENT] Error embedding chunks for content %s: %r",
            content_id[:8],
            e,
        )
        raise

    # Prepare vectors for Pinecone
    cfg = get_config()
    store = get_vector_store()

    vectors: List[Dict[str, Any]] = []
    embedding_ids: List[str] = []

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{content_id}_chunk_{order:03d}"
        embedding_ids.append(chunk_id)

        snippet = textwrap.shorten(chunk, width=300)

        meta: Dict[str, Any] = {
            **base_meta,
            "chunk_id": chunk_id,
            "chunk_index": order,
            "snippet": snippet,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(meta),
            }
        )

    if not vectors:
        logger.warning(
            "[INGEST][POSTGRES-CONTENT] No vectors to upsert for content %s",
            content_id[:8],
        )
        return []

    logger.info(
        "[INGEST][POSTGRES-CONTENT] Upserting %d vectors into Pinecone for content %s",
        len(vectors),
        content_id[:8],
    )

    # Batch upserts
    batch_size = 100
    try:
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(vectors) + batch_size - 1) // batch_size

            store.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

            if batch_num % 10 == 0 or batch_num == total_batches:
                logger.info(
                    "[INGEST][POSTGRES-CONTENT] Completed batch %d/%d for content %s",
                    batch_num,
                    total_batches,
                    content_id[:8],
                )
    except Exception as e:
        logger.error(
            "[INGEST][POSTGRES-CONTENT] Error upserting vectors for content %s: %r",
            content_id[:8],
            e,
        )
        raise

    # Extract features and detect signatures (optional, non-blocking)
    try:
        feature_names = extract_features_from_text(content)
        if feature_names:
            logger.info(
                "[INGEST][POSTGRES-CONTENT] Extracted %d feature(s) from content %s",
                len(feature_names),
                content_id[:8],
            )
            # Note: Features could be linked to datasets if needed
            # For now, just log them
    except Exception as e:
        logger.debug(
            "[INGEST][POSTGRES-CONTENT] Error extracting features: %r",
            e,
        )

    # Detect signatures (optional, requires Notion page ID - skip if not available)
    try:
        {
            "diseases": base_meta.get("diseases", []),
            "matrix": base_meta.get("matrix", []),
            "model_systems": base_meta.get("model_systems", []),
        }

        # Signature detection requires Notion page ID - skip for Postgres-only mode
        logger.debug(
            "[INGEST][POSTGRES-CONTENT] Skipping signature detection (Postgres-only mode)",
        )
    except Exception as e:
        logger.debug(
            "[INGEST][POSTGRES-CONTENT] Error in signature detection (skipped): %r",
            e,
        )

    logger.info(
        "[INGEST][POSTGRES-CONTENT] Direct ingestion complete for content %s",
        content_id[:8],
    )

    return embedding_ids


def ingest_email_content(
    email_content: str,
    title: str,
    from_sender: Optional[str] = None,
    email_id: Optional[str] = None,
    tags: Optional[List[str]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    force: bool = False,
) -> List[str]:
    """
    Ingest email content directly to Pinecone (Postgres-only, no Notion).

    Idempotent: Only ingests if email is new or changed.

    Args:
        email_content: Full email text content
        title: Email subject/title
        from_sender: Sender email address
        email_id: Optional unique email identifier (required for idempotency)
        tags: Optional list of tags
        metadata: Optional additional metadata
        force: If True, re-ingest even if already ingested

    Returns:
        List of embedding/chunk IDs created (empty if skipped)
    """
    # Build email metadata
    email_meta: Dict[str, Any] = {
        "doc_type": "Email",
        "doc_source": "Email",
        "from": from_sender[:200] if from_sender else "",
        "tags": tags or [],
        **(metadata or {}),
    }

    # Add email header to content
    header = f"Title: {title}\n"
    if from_sender:
        header += f"From: {from_sender}\n"
    header += "\n"

    full_content = header + email_content

    return ingest_content_direct_to_pinecone(
        content=full_content,
        source_type="email",
        title=title,
        metadata=email_meta,
        content_id=email_id,
        force=force,
    )


def ingest_literature_content(
    literature_content: str,
    title: str,
    authors: Optional[List[str]] = None,
    doi: Optional[str] = None,
    zotero_key: Optional[str] = None,
    metadata: Optional[Dict[str, Any]] = None,
    force: bool = False,
) -> List[str]:
    """
    Ingest literature/Zotero content directly to Pinecone (Postgres-only, no Notion).

    Idempotent: Only ingests if literature is new or changed.

    Args:
        literature_content: Full literature text content
        title: Publication title
        authors: Optional list of authors
        doi: Optional DOI
        zotero_key: Optional Zotero item key (used as content_id for idempotency)
        metadata: Optional additional metadata
        force: If True, re-ingest even if already ingested

    Returns:
        List of embedding/chunk IDs created (empty if skipped)
    """
    # Build literature metadata
    lit_meta: Dict[str, Any] = {
        "doc_type": "Literature",
        "doc_source": "Zotero",
        "authors": authors or [],
        "doi": doi or "",
        "zotero_key": zotero_key or "",
        **(metadata or {}),
    }

    # Add literature header to content
    header = f"Title: {title}\n"
    if authors:
        header += f"Authors: {', '.join(authors[:5])}\n"  # Limit to first 5 authors
    if doi:
        header += f"DOI: {doi}\n"
    header += "\n"

    full_content = header + literature_content

    return ingest_content_direct_to_pinecone(
        content=full_content,
        source_type="literature",
        title=title,
        metadata=lit_meta,
        content_id=zotero_key or None,
        force=force,
    )

