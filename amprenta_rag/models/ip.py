"""IP and patent tracking models."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import TYPE_CHECKING, List, Optional

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.models.auth import Company, User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class InventionDisclosure(Base):
    """Invention disclosure for IP tracking."""

    __tablename__ = "invention_disclosures"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(500), nullable=False)
    description = Column(Text, nullable=False)
    technical_field = Column(String(200), nullable=True)
    status = Column(String(50), nullable=False, default="draft")  # draft/submitted/under_review/filed/granted/rejected
    priority_date = Column(DateTime, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    company_id = Column(UUID(as_uuid=True), ForeignKey("companies.id"), nullable=False, index=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    __table_args__ = (
        Index("ix_invention_disclosures_company_status", "company_id", "status"),
        Index("ix_invention_disclosures_priority_date", "priority_date"),
    )

    created_by: Mapped["User"] = relationship(foreign_keys=[created_by_id])
    company: Mapped["Company"] = relationship(foreign_keys=[company_id])
    inventors: Mapped[List["DisclosureInventor"]] = relationship(back_populates="disclosure", cascade="all, delete-orphan")
    patents: Mapped[List["PatentApplication"]] = relationship(back_populates="disclosure", cascade="all, delete-orphan")
    links: Mapped[List["IPLink"]] = relationship(back_populates="disclosure", cascade="all, delete-orphan")


class DisclosureInventor(Base):
    """Inventor relationship to invention disclosure."""

    __tablename__ = "disclosure_inventors"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    disclosure_id = Column(UUID(as_uuid=True), ForeignKey("invention_disclosures.id"), nullable=False, index=True)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    contribution_percentage = Column(Integer, nullable=True)
    is_primary = Column(Boolean, default=False, nullable=False)

    disclosure: Mapped["InventionDisclosure"] = relationship(back_populates="inventors")
    user: Mapped["User"] = relationship(foreign_keys=[user_id])


class PatentApplication(Base):
    """Patent application linked to invention disclosure."""

    __tablename__ = "patent_applications"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    disclosure_id = Column(UUID(as_uuid=True), ForeignKey("invention_disclosures.id"), nullable=False, index=True)
    application_number = Column(String(100), nullable=False)
    filing_date = Column(DateTime, nullable=False)
    jurisdiction = Column(String(50), nullable=False)  # US, EP, PCT, etc.
    status = Column(String(50), nullable=False, default="pending", index=True)
    grant_date = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    disclosure: Mapped["InventionDisclosure"] = relationship(back_populates="patents")
    claims: Mapped[List["PatentClaim"]] = relationship(back_populates="patent", cascade="all, delete-orphan")


class PatentClaim(Base):
    """Individual patent claim."""

    __tablename__ = "patent_claims"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    patent_id = Column(UUID(as_uuid=True), ForeignKey("patent_applications.id"), nullable=False, index=True)
    claim_number = Column(Integer, nullable=False)
    claim_text = Column(Text, nullable=False)
    claim_type = Column(String(20), nullable=False)  # independent/dependent
    parent_claim_id = Column(UUID(as_uuid=True), ForeignKey("patent_claims.id"), nullable=True)

    patent: Mapped["PatentApplication"] = relationship(back_populates="claims")
    parent_claim: Mapped[Optional["PatentClaim"]] = relationship(remote_side=[id])


class IPLink(Base):
    """Link between entities and invention disclosures."""

    __tablename__ = "ip_links"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    entity_type = Column(String(50), nullable=False, index=True)
    entity_id = Column(UUID(as_uuid=True), nullable=False, index=True)
    disclosure_id = Column(UUID(as_uuid=True), ForeignKey("invention_disclosures.id"), nullable=False, index=True)
    link_type = Column(String(50), nullable=False)
    notes = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    __table_args__ = (
        UniqueConstraint("entity_type", "entity_id", "disclosure_id", "link_type", name="uq_ip_links"),
        Index("ix_ip_links_entity", "entity_type", "entity_id"),
    )

    disclosure: Mapped["InventionDisclosure"] = relationship(back_populates="links")


__all__ = [
    "InventionDisclosure",
    "DisclosureInventor",
    "PatentApplication",
    "PatentClaim",
    "IPLink",
]

