"""
Embedding pipeline for scientific papers.

Creates section-aware RAG chunks from parsed paper content and stores
them with embeddings for semantic search.
"""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.ingestion.papers.jats_parser import PaperContent
from amprenta_rag.ingestion.postgres_rag_chunk import create_rag_chunk_in_postgres
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_paper_chunks(
    paper_content: PaperContent,
    literature_id: UUID,
    pmid: str,
    paper_title: str,
    db: Session,
) -> List[UUID]:
    """
    Create RAG chunks from paper sections.

    Creates one chunk per section with section title as header.
    Chunks are linked to the Literature record for retrieval.

    Args:
        paper_content: Parsed paper content with sections
        literature_id: UUID of Literature record
        pmid: PubMed ID
        paper_title: Paper title
        db: Database session

    Returns:
        List of created chunk UUIDs

    Example:
        >>> from amprenta_rag.ingestion.papers import parse_jats_xml
        >>> jats_xml = "<article>...</article>"
        >>> content = parse_jats_xml(jats_xml)
        >>> chunk_ids = create_paper_chunks(content, lit_id, "12345", "Title", db)
        >>> len(chunk_ids) == len(content.sections)
        True
    """
    chunk_ids = []

    for idx, section in enumerate(paper_content.sections):
        # Create chunk ID from PMID and section index
        chunk_id = f"pmid_{pmid}_section_{idx:03d}"

        # Combine section title and content
        chunk_text = f"{section.title}\n\n{section.content}"

        # Create snippet (first 500 chars)
        snippet = chunk_text[:500] if len(chunk_text) > 500 else chunk_text

        # Metadata for this chunk
        chunk_metadata = {
            "section_title": section.title,
            "section_order": section.order,
            "pmid": pmid,
            "paper_title": paper_title,
        }

        try:
            chunk_uuid = create_rag_chunk_in_postgres(
                chunk_id=chunk_id,
                chunk_text=chunk_text,
                source_type="Literature",
                source_id=literature_id,
                source_name=paper_title,
                chunk_index=idx,
                snippet=snippet,
                chunk_metadata=chunk_metadata,
                literature_id=literature_id,
                db=db,
            )
            chunk_ids.append(chunk_uuid)
            logger.debug(
                "[PAPER_EMBED] Created chunk %d/%d for PMID %s: %s",
                idx + 1,
                len(paper_content.sections),
                pmid,
                section.title,
            )
        except Exception as e:
            logger.error(
                "[PAPER_EMBED] Failed to create chunk for section '%s' (PMID %s): %r",
                section.title,
                pmid,
                e,
            )
            # Continue with other sections
            continue

    logger.info(
        "[PAPER_EMBED] Created %d chunks for PMID %s (%s)",
        len(chunk_ids),
        pmid,
        paper_title,
    )
    return chunk_ids


def embed_paper_sections(
    paper_content: PaperContent,
    literature_id: UUID,
    pmid: str,
    paper_title: str,
    db: Optional[Session] = None,
) -> List[UUID]:
    """
    Create and embed paper sections as RAG chunks.

    Convenience function that wraps create_paper_chunks with optional
    session management.

    Args:
        paper_content: Parsed paper content
        literature_id: UUID of Literature record
        pmid: PubMed ID
        paper_title: Paper title
        db: Optional database session

    Returns:
        List of created chunk UUIDs
    """
    if db is None:
        from amprenta_rag.database.base import get_db

        db_gen = get_db()
        session = next(db_gen)
        try:
            return create_paper_chunks(
                paper_content, literature_id, pmid, paper_title, session
            )
        finally:
            db_gen.close()
    else:
        return create_paper_chunks(paper_content, literature_id, pmid, paper_title, db)

