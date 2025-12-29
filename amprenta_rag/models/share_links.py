"""Share link models for Voila dashboard sharing."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Index, Integer, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.models.auth import Company, User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class ShareLink(Base):
    """Shareable link for Voila dashboards."""

    __tablename__ = "share_links"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    token = Column(String(64), unique=True, nullable=False, index=True)
    dashboard_path = Column(String(500), nullable=False)
    context_json = Column(Text, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False)
    company_id = Column(UUID(as_uuid=True), ForeignKey("companies.id"), nullable=False, index=True)
    expires_at = Column(DateTime, nullable=False)
    max_views = Column(Integer, nullable=True)
    view_count = Column(Integer, default=0, nullable=False)
    is_active = Column(Boolean, default=True, nullable=False)
    permissions = Column(String(50), default="view", nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    last_accessed_at = Column(DateTime, nullable=True)

    __table_args__ = (
        Index("ix_share_links_company_active", "company_id", "is_active"),
        Index("ix_share_links_expires", "expires_at"),
    )

    created_by: Mapped["User"] = relationship(foreign_keys=[created_by_id])
    company: Mapped["Company"] = relationship(foreign_keys=[company_id])


__all__ = ["ShareLink"]

