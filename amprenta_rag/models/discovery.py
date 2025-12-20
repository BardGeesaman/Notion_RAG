"""Discovery and harvesting models."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import Optional, TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, ForeignKey, Integer, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import Experiment, User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class DiscoveryJob(Base):
    """Track automated repository discovery jobs."""

    __tablename__ = "discovery_jobs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    repository = Column(String(50), nullable=False)
    query = Column(String(500), nullable=False)
    status = Column(String(50), default="pending")
    studies_found = Column(Integer, default=0)
    studies_imported = Column(Integer, default=0)
    error_message = Column(Text, nullable=True)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)

    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class DiscoveredStudy(Base):
    """Studies found during discovery jobs."""

    __tablename__ = "discovered_studies"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    job_id = Column(UUID(as_uuid=True), ForeignKey("discovery_jobs.id"), nullable=False)
    study_id = Column(String(100), nullable=False)
    repository = Column(String(50), nullable=False)
    title = Column(String(500), nullable=True)
    description = Column(Text, nullable=True)
    organism = Column(String(200), nullable=True)
    omics_type = Column(String(100), nullable=True)
    status = Column(String(50), default="new")
    imported_experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True)
    discovered_at = Column(DateTime, default=datetime.utcnow)

    job: Mapped["DiscoveryJob"] = relationship("DiscoveryJob", backref="discovered_studies")
    imported_experiment: Mapped[Optional["Experiment"]] = relationship("Experiment")


class HarvestSchedule(Base):
    """Scheduled automated harvesting from repositories."""

    __tablename__ = "harvest_schedules"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    repository = Column(String(50), nullable=False)
    query = Column(String(500), nullable=False)
    interval_hours = Column(Integer, default=24)
    is_active = Column(Boolean, default=True)
    last_run = Column(DateTime(timezone=True), nullable=True)
    next_run = Column(DateTime(timezone=True), nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


__all__ = ["DiscoveryJob", "DiscoveredStudy", "HarvestSchedule"]

