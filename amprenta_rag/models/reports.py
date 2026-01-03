"""Report Builder models."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import TYPE_CHECKING, Optional
from uuid import uuid4

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, JSON, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import Program, User


def generate_uuid():
    return uuid4()


class ReportTemplate(Base):
    """User-created report template with configurable sections."""
    
    __tablename__ = "report_templates"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    
    # Sections stored as JSON array: [{type, config, order}]
    # Validated by SECTION_SCHEMAS before save
    sections = Column(JSON, nullable=False, default=list)
    
    # Sharing
    is_public = Column(Boolean, default=False, nullable=False)
    
    # Scoping
    program_id = Column(UUID(as_uuid=True), ForeignKey("programs.id"), nullable=True)
    
    # Audit
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime(timezone=True), onupdate=lambda: datetime.now(timezone.utc))
    
    # Relationships
    program: Mapped[Optional["Program"]] = relationship("Program", foreign_keys=[program_id])
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    
    def __repr__(self) -> str:
        return f"<ReportTemplate(id={self.id}, name='{self.name}')>"
