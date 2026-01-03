"""Data Catalog models for metadata management."""

from datetime import datetime
from typing import List, Optional
from uuid import uuid4

from sqlalchemy import (
    Boolean, DateTime, ForeignKey, Index, Integer, String, Text,
    UniqueConstraint, func
)
from sqlalchemy.dialects.postgresql import JSONB, UUID as PostgresUUID
from sqlalchemy.orm import Mapped, mapped_column, relationship

from amprenta_rag.database.base import Base


class CatalogEntry(Base):
    """Registered entity type (table) in the catalog."""
    __tablename__ = "catalog_entries"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    entity_type: Mapped[str] = mapped_column(String(100), unique=True)
    table_name: Mapped[str] = mapped_column(String(100))
    display_name: Mapped[str] = mapped_column(String(200))
    description: Mapped[Optional[str]] = mapped_column(Text)
    category: Mapped[str] = mapped_column(String(50))  # Core, Chemistry, Omics, Admin
    row_count: Mapped[Optional[int]] = mapped_column(Integer)
    last_refreshed: Mapped[Optional[datetime]] = mapped_column(DateTime)
    is_active: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    columns: Mapped[List["ColumnMetadata"]] = relationship(back_populates="catalog_entry", cascade="all, delete-orphan")


class ColumnMetadata(Base):
    """Column-level metadata for catalog entries."""
    __tablename__ = "column_metadata"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    catalog_entry_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("catalog_entries.id", ondelete="CASCADE"))
    column_name: Mapped[str] = mapped_column(String(100))
    display_name: Mapped[Optional[str]] = mapped_column(String(200))
    data_type: Mapped[str] = mapped_column(String(50))
    description: Mapped[Optional[str]] = mapped_column(Text)
    is_nullable: Mapped[bool] = mapped_column(Boolean, default=True)
    is_primary_key: Mapped[bool] = mapped_column(Boolean, default=False)
    is_foreign_key: Mapped[bool] = mapped_column(Boolean, default=False)
    foreign_key_target: Mapped[Optional[str]] = mapped_column(String(200))
    example_values: Mapped[Optional[list]] = mapped_column(JSONB)  # P1 FIX: JSONB
    glossary_term_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("glossary_terms.id", ondelete="SET NULL"))
    
    catalog_entry: Mapped["CatalogEntry"] = relationship(back_populates="columns")
    glossary_term: Mapped[Optional["GlossaryTerm"]] = relationship()
    
    __table_args__ = (
        UniqueConstraint('catalog_entry_id', 'column_name', name='uix_column_entry_name'),
    )


class GlossaryTerm(Base):
    """Business glossary terms with definitions."""
    __tablename__ = "glossary_terms"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    term: Mapped[str] = mapped_column(String(200), unique=True)
    definition: Mapped[str] = mapped_column(Text)
    category: Mapped[Optional[str]] = mapped_column(String(100))
    synonyms: Mapped[Optional[list]] = mapped_column(JSONB)  # P1 FIX: JSONB
    related_terms: Mapped[Optional[list]] = mapped_column(JSONB)  # P1 FIX: JSONB
    source: Mapped[Optional[str]] = mapped_column(String(200))
    created_by_id: Mapped[Optional[uuid4]] = mapped_column(PostgresUUID(as_uuid=True), ForeignKey("users.id", ondelete="SET NULL"))
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=func.now(), onupdate=func.now())


class DataLineageEdge(Base):
    """Data flow between entities."""
    __tablename__ = "data_lineage_edges"
    
    id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True), primary_key=True, default=uuid4)
    source_type: Mapped[str] = mapped_column(String(100))
    source_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True))
    target_type: Mapped[str] = mapped_column(String(100))
    target_id: Mapped[uuid4] = mapped_column(PostgresUUID(as_uuid=True))
    relationship_type: Mapped[str] = mapped_column(String(50))  # derived_from, input_to, generated_by
    transformation: Mapped[Optional[str]] = mapped_column(String(200))
    edge_metadata: Mapped[Optional[dict]] = mapped_column(JSONB)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=func.now())
    
    __table_args__ = (
        UniqueConstraint('source_type', 'source_id', 'target_type', 'target_id', 
                         name='uix_lineage_source_target'),  # P1 FIX
        Index("ix_lineage_source", "source_type", "source_id"),
        Index("ix_lineage_target", "target_type", "target_id"),
    )


__all__ = ["CatalogEntry", "ColumnMetadata", "GlossaryTerm", "DataLineageEdge"]
