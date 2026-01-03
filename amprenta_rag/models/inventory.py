"""Compound inventory models for physical sample tracking."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import Optional, TYPE_CHECKING

from sqlalchemy import Column, DateTime, Float, ForeignKey, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import User
    from amprenta_rag.models.sample import Sample, StorageLocation


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class CompoundPlate(Base):
    """Compound storage plate (mother/daughter plates)."""

    __tablename__ = "compound_plates"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    barcode = Column(String(200), nullable=False, unique=True, index=True)
    plate_format = Column(String(20), nullable=False)  # 96, 384, 1536
    plate_type = Column(String(50), nullable=False)  # mother, daughter, screening
    storage_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    status = Column(String(20), nullable=False, default="active")  # active, archived, empty
    created_at = Column(DateTime, default=datetime.utcnow)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    storage_location: Mapped[Optional["StorageLocation"]] = relationship("StorageLocation", backref="compound_plates")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class CompoundRequest(Base):
    """Request for compound samples from inventory."""

    __tablename__ = "compound_requests"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id"), nullable=True, index=True)
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)
    requester_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=False, index=True)
    
    requested_quantity = Column(Float, nullable=False)
    quantity_unit = Column(String(20), nullable=False, default="µL")
    purpose = Column(String(100), nullable=True)  # screening, assay, synthesis
    priority = Column(String(20), nullable=False, default="normal")  # low, normal, high, urgent
    status = Column(String(20), nullable=False, default="requested")  # requested, approved, fulfilled, cancelled, rejected
    
    requested_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    approved_at = Column(DateTime, nullable=True)
    approved_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    fulfilled_at = Column(DateTime, nullable=True)
    fulfilled_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    rejection_reason = Column(Text, nullable=True)
    notes = Column(Text, nullable=True)

    sample: Mapped[Optional["Sample"]] = relationship("Sample", backref="requests")
    compound: Mapped[Optional["Compound"]] = relationship("Compound", backref="requests")
    requester: Mapped["User"] = relationship("User", foreign_keys=[requester_id], backref="compound_requests")
    approved_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[approved_by_id])
    fulfilled_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[fulfilled_by_id])


__all__ = ["CompoundPlate", "CompoundRequest"]
