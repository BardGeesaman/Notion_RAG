"""Electronic lab notebook and protocol models."""

from __future__ import annotations

import uuid
from datetime import datetime

from sqlalchemy import ARRAY, Boolean, Column, DateTime, ForeignKey, Integer, JSON, String, Text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship

from amprenta_rag.database.base import Base


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class LabNotebookEntry(Base):
    """Electronic Lab Notebook entry for tracking experiments, observations, and notes."""

    __tablename__ = "lab_notebook_entries"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(500), nullable=False)
    content = Column(Text, nullable=False)
    entry_type = Column(String(50), nullable=True)
    tags = Column(ARRAY(String), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False, index=True)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    attachments = Column(JSON, nullable=True)
    entry_metadata = Column(JSON, nullable=True)
    notion_page_id = Column(String(36), nullable=True, unique=True, index=True)


class LabNotebookEntryAssociation(Base):
    """Association model for linking lab notebook entries to various entities."""

    __tablename__ = "lab_notebook_entry_assoc"

    entry_id = Column(UUID(as_uuid=True), ForeignKey("lab_notebook_entries.id"), primary_key=True)
    entity_type = Column(String(50), primary_key=True)
    entity_id = Column(UUID(as_uuid=True), primary_key=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    entry = relationship("LabNotebookEntry", backref="linked_entities")


class Protocol(Base):
    """Reusable experimental protocol template."""

    __tablename__ = "protocols"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    version = Column(Integer, default=1, nullable=False)
    category = Column(String(100), nullable=True)
    description = Column(Text, nullable=True)
    steps = Column(JSON, nullable=True)
    materials = Column(JSON, nullable=True)
    parameters = Column(JSON, nullable=True)
    is_template = Column(Boolean, default=True)
    is_active = Column(Boolean, default=True)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("protocols.id"), nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    parent = relationship("Protocol", remote_side=[id], backref="versions")
    created_by = relationship("User", foreign_keys=[created_by_id])


class ExperimentProtocol(Base):
    """Link between experiments and protocols used."""

    __tablename__ = "experiment_protocols"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=False)
    protocol_id = Column(UUID(as_uuid=True), ForeignKey("protocols.id"), nullable=False)
    executed_at = Column(DateTime, nullable=True)
    notes = Column(Text, nullable=True)
    deviations = Column(JSON, nullable=True)

    experiment = relationship("Experiment", backref="protocol_links")
    protocol = relationship("Protocol", backref="experiment_links")


__all__ = ["LabNotebookEntry", "LabNotebookEntryAssociation", "Protocol", "ExperimentProtocol"]

