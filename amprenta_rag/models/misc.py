"""Miscellaneous domain models (templates, scheduling, retention, ontology, variants, costs)."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import Optional, TYPE_CHECKING

from sqlalchemy import ARRAY, Boolean, Column, DateTime, Float, ForeignKey, Integer, JSON, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.database.models import Experiment, Project, User


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class ExperimentTemplate(Base):
    """Reusable experiment templates for quick creation."""

    __tablename__ = "experiment_templates"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    design_type = Column(String(50), nullable=True)
    organism = Column(ARRAY(String), nullable=True)
    sample_groups = Column(JSON, nullable=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class ScheduledEvent(Base):
    """Scheduled events for experiments and resources."""

    __tablename__ = "scheduled_events"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    title = Column(String(500), nullable=False)
    event_type = Column(String(100), nullable=False, index=True)  # e.g., "experiment", "equipment", "meeting", "maintenance"
    resource_name = Column(String(200), nullable=False, index=True)  # Equipment name, room, etc.
    start_time = Column(DateTime, nullable=False, index=True)
    end_time = Column(DateTime, nullable=False, index=True)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True, index=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    notes = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    experiment: Mapped[Optional["Experiment"]] = relationship("Experiment", backref="scheduled_events")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class RetentionPolicy(Base):
    """Data retention policies for entities."""

    __tablename__ = "retention_policies"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(200), nullable=False)
    entity_type = Column(String(50), nullable=False, index=True)  # "experiment", "dataset", "compound", etc.
    retention_days = Column(Integer, nullable=False)  # Number of days to retain before action
    action = Column(String(50), nullable=False)  # "archive", "delete", "notify"
    is_active = Column(Boolean, default=True, nullable=False, index=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


class OntologyTerm(Base):
    """Ontology terms for controlled vocabularies."""

    __tablename__ = "ontology_terms"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    vocabulary = Column(String(100), nullable=False, index=True)  # e.g., "disease", "tissue", "cell_type"
    term = Column(String(500), nullable=False, index=True)
    description = Column(Text, nullable=True)
    parent_id = Column(UUID(as_uuid=True), ForeignKey("ontology_terms.id"), nullable=True, index=True)
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    parent: Mapped[Optional["OntologyTerm"]] = relationship(
        "OntologyTerm", remote_side=[id], backref="children"
    )


class GeneticVariant(Base):
    """Genetic variants tracked in experiments and cell lines."""

    __tablename__ = "genetic_variants"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    gene = Column(String, nullable=False, index=True)
    variant = Column(String, nullable=False)
    zygosity = Column(String)  # homozygous, heterozygous, hemizygous
    organism = Column(String, nullable=False)  # cell line or organism name
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True)
    notes = Column(Text)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    
    # VCF-specific fields (nullable for backward compatibility)
    chromosome = Column(String(50), nullable=True, index=True)
    position = Column(Integer, nullable=True)
    reference_allele = Column(String(500), nullable=True)
    alternate_allele = Column(String(500), nullable=True)
    quality = Column(Float, nullable=True)
    vcf_filter = Column(String(100), nullable=True)

    experiment: Mapped[Optional["Experiment"]] = relationship("Experiment", backref="variants")


class AlignmentFile(Base):
    """BAM/CRAM alignment file metadata and statistics."""

    __tablename__ = "alignment_files"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    filename = Column(String(500), nullable=False)
    file_path = Column(String(1000), nullable=False)
    index_file_path = Column(String(1000), nullable=True)
    file_format = Column(String(10), nullable=False)  # "BAM" or "CRAM"
    has_index = Column(Boolean, default=False, nullable=False, index=True)
    
    # Header metadata
    reference_genome = Column(String(100), nullable=True)
    num_references = Column(Integer, nullable=True)
    read_groups = Column(JSON, nullable=True)
    
    # Stats
    total_reads = Column(Integer, nullable=True)  # Use Integer, BigInteger if needed
    mapped_reads = Column(Integer, nullable=True)
    unmapped_reads = Column(Integer, nullable=True)
    duplicate_rate = Column(Float, nullable=True)
    mean_coverage = Column(Float, nullable=True)
    
    # Relationships
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True, index=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

    experiment: Mapped[Optional["Experiment"]] = relationship("Experiment", backref="alignment_files")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


class CostEntry(Base):
    """Cost tracking entries for projects and experiments."""

    __tablename__ = "cost_entries"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    project_id = Column(UUID(as_uuid=True), ForeignKey("projects.id"), nullable=True, index=True)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("experiments.id"), nullable=True, index=True)
    category = Column(String(100), nullable=False, index=True)  # e.g., "equipment", "reagents", "personnel", "overhead"
    description = Column(Text, nullable=False)
    amount = Column(Float, nullable=False)
    currency = Column(String(10), nullable=False, default="USD")
    entry_date = Column(DateTime, nullable=False, index=True)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    project: Mapped[Optional["Project"]] = relationship("Project", backref="cost_entries")
    experiment: Mapped[Optional["Experiment"]] = relationship("Experiment", backref="cost_entries")
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])


__all__ = [
    "ExperimentTemplate",
    "ScheduledEvent",
    "RetentionPolicy",
    "OntologyTerm",
    "GeneticVariant",
    "AlignmentFile",
    "CostEntry",
]

