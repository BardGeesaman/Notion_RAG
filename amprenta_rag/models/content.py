"""Content models: literature, emails, RAG chunks, critiques."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Integer, JSON, String, Text, func
from sqlalchemy.dialects.postgresql import ARRAY, UUID, TSVECTOR
from sqlalchemy.orm import Mapped, relationship

try:
    from pgvector.sqlalchemy import Vector
except ImportError:
    # Fallback for environments without pgvector installed
    Vector = None  # type: ignore

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class Literature(Base):
    """Literature item (from Zotero or other sources)."""

    __tablename__ = "literature"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(1000), nullable=False, index=True)
    source_type = Column(String(100), nullable=True)
    abstract = Column(Text, nullable=True)
    doi = Column(String(500), nullable=True, index=True)
    url = Column(String(1000), nullable=True)
    date = Column(String(100), nullable=True)
    journal = Column(String(500), nullable=True)
    year = Column(Integer, nullable=True)
    tags = Column(ARRAY(String), nullable=True)
    zotero_item_key = Column(String(50), nullable=True, unique=True, index=True)
    zotero_item_type = Column(String(100), nullable=True)
    embedding_status = Column(String(50), nullable=True, default="Not Embedded")
    last_ingested_at = Column(DateTime, nullable=True)
    semantic_metadata = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)
    external_ids = Column(JSON, nullable=True)

    # Paper source tracking
    source = Column(String(50), default="zotero", nullable=False)
    
    # PubMed/PMC identifiers
    pmid = Column(String(20), nullable=True, index=True)
    pmc_id = Column(String(20), nullable=True)
    
    # MeSH terms for faceted search
    mesh_terms = Column(ARRAY(String), default=list)
    
    # Full text availability flag
    full_text_available = Column(Boolean, default=False)
    
    # Semantic Scholar / OpenAlex integration
    semantic_scholar_id = Column(String(50), nullable=True, index=True)
    openalex_id = Column(String(50), nullable=True, index=True)
    tldr_summary = Column(Text, nullable=True)  # AI-generated summary from S2
    citation_count = Column(Integer, nullable=True)
    influential_citation_count = Column(Integer, nullable=True)
    open_access_status = Column(String(50), nullable=True)

    chunks: Mapped[List["RAGChunk"]] = relationship(
        back_populates="literature", cascade="all, delete-orphan"
    )
    
    # Citation relationships
    citations_made: Mapped[List["PaperCitation"]] = relationship(
        foreign_keys="[PaperCitation.citing_paper_id]",
        back_populates="citing_paper",
        cascade="all, delete-orphan"
    )
    
    citations_received: Mapped[List["PaperCitation"]] = relationship(
        foreign_keys="[PaperCitation.cited_paper_id]",
        back_populates="cited_paper",
        cascade="all, delete-orphan"
    )


class Email(Base):
    """Email or Note item."""

    __tablename__ = "emails"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(1000), nullable=False, index=True)
    from_sender = Column(String(500), nullable=True)
    item_type = Column(String(100), nullable=True)  # Email, Note, etc.
    tags = Column(ARRAY(String), nullable=True)
    content = Column(Text, nullable=True)
    semantic_metadata = Column(JSON, nullable=True)
    embedding_status = Column(String(50), nullable=True, default="Not Embedded")
    last_ingested_at = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)
    external_ids = Column(JSON, nullable=True)

    chunks: Mapped[List["RAGChunk"]] = relationship(
        back_populates="email", cascade="all, delete-orphan"
    )


class RAGChunk(Base):
    """RAG chunk from literature, datasets, or other sources."""

    __tablename__ = "rag_chunks"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    chunk_id = Column(String(500), nullable=False, unique=True, index=True)
    chunk_text = Column(Text, nullable=False)
    chunk_index = Column(Integer, nullable=False, default=0)
    snippet = Column(String(500), nullable=True)
    source_type = Column(String(100), nullable=False)
    source_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    source_name = Column(String(500), nullable=True)
    zotero_item_key = Column(String(50), nullable=True, index=True)
    note_key = Column(String(50), nullable=True)
    attachment_key = Column(String(50), nullable=True)
    note_hash = Column(String(32), nullable=True)
    attachment_hash = Column(String(32), nullable=True)
    chunk_metadata = Column(JSON, nullable=True)
    search_vector = Column(TSVECTOR, nullable=True)
    embedding = Column(Vector(3072) if Vector is not None else Text, nullable=True)
    embedding_model = Column(String(100), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)

    literature_id = Column(UUID(as_uuid=True), ForeignKey("literature.id"), nullable=True, index=True)
    literature: Mapped[Optional["Literature"]] = relationship(back_populates="chunks")

    email_id = Column(UUID(as_uuid=True), ForeignKey("emails.id"), nullable=True, index=True)
    email: Mapped[Optional["Email"]] = relationship(back_populates="chunks")

    extraction_job_id = Column(
        UUID(as_uuid=True), ForeignKey("extraction_jobs.id"), nullable=True, index=True
    )


class PaperCitation(Base):
    """Citation relationship between papers."""

    __tablename__ = "paper_citations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    
    # Citing paper (the paper that cites another)
    citing_paper_id = Column(UUID(as_uuid=True), ForeignKey("literature.id"), nullable=False, index=True)
    citing_paper: Mapped["Literature"] = relationship(
        foreign_keys=[citing_paper_id],
        back_populates="citations_made"
    )
    
    # Cited paper (the paper being cited - may not be in our system yet)
    cited_paper_id = Column(UUID(as_uuid=True), ForeignKey("literature.id"), nullable=True, index=True)
    cited_paper: Mapped[Optional["Literature"]] = relationship(
        foreign_keys=[cited_paper_id],
        back_populates="citations_received"
    )
    
    # Citation metadata for papers not yet in system
    cited_doi = Column(String(500), nullable=True)
    cited_title = Column(String(1000), nullable=True)
    
    # Semantic Scholar specific fields
    citation_context = Column(Text, nullable=True)  # Context snippet from citing paper
    is_influential = Column(Boolean, default=False)  # S2's influential citation flag
    
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


class LiteratureCritique(Base):
    """Critical analysis of scientific literature."""

    __tablename__ = "literature_critiques"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    source_type = Column(String(50), nullable=False)
    source_id = Column(UUID(as_uuid=True), nullable=True)
    source_text = Column(Text, nullable=False)
    critique_type = Column(String(50), nullable=False)
    content = Column(JSON, nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by: Mapped[Optional["User"]] = relationship(foreign_keys=[created_by_id])


__all__ = ["Literature", "Email", "RAGChunk", "PaperCitation", "LiteratureCritique"]

