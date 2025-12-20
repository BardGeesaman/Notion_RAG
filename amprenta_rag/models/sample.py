"""Sample inventory models."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import Experiment, User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class StorageLocation(Base):
    """Physical storage hierarchy for samples (freezer/shelf/box/rack)."""

    __tablename__ = "storage_locations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    location_type = Column(String(50), nullable=True)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    temperature = Column(String(50), nullable=True)
    capacity = Column(Integer, nullable=True)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    parent: Mapped[Optional["StorageLocation"]] = relationship(
        "StorageLocation", remote_side=[id], backref="children"
    )


class Sample(Base):
    """Biological sample with storage and lineage tracking."""

    __tablename__ = "samples"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    sample_type = Column(String(100), nullable=True)
    barcode = Column(String(255), unique=True, nullable=True, index=True)
    storage_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    position = Column(String(100), nullable=True)
    parent_sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id"), nullable=True)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True)
    quantity = Column(Float, nullable=True)
    unit = Column(String(50), nullable=True)
    status = Column(String(50), default="available")
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    notes = Column(Text, nullable=True)

    storage_location: Mapped[Optional["StorageLocation"]] = relationship("StorageLocation", backref="samples")
    parent_sample: Mapped[Optional["Sample"]] = relationship(
        "Sample", remote_side=[id], backref="child_samples"
    )
    experiment: Mapped[Optional["Experiment"]] = relationship("Experiment", backref="samples")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class SampleTransfer(Base):
    """Track movement of samples between storage locations."""

    __tablename__ = "sample_transfers"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id"), nullable=False)
    from_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    to_location_id = Column(UUID(as_uuid=True), ForeignKey("storage_locations.id"), nullable=True)
    transferred_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    transferred_at = Column(DateTime, default=datetime.utcnow)
    notes = Column(Text, nullable=True)

    sample: Mapped["Sample"] = relationship("Sample", backref="transfers")
    from_location: Mapped[Optional["StorageLocation"]] = relationship(
        "StorageLocation", foreign_keys=[from_location_id], backref="outgoing_transfers"
    )
    to_location: Mapped[Optional["StorageLocation"]] = relationship(
        "StorageLocation", foreign_keys=[to_location_id], backref="incoming_transfers"
    )
    transferred_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[transferred_by_id])


__all__ = ["StorageLocation", "Sample", "SampleTransfer"]

