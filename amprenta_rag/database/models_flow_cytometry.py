"""
Flow Cytometry / FACS database models for the multi-omics platform.

These models store flow cytometry dataset metadata, parameters, gates, and population statistics.
Event-level data is stored in Parquet files for performance.
"""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from enum import Enum as PyEnum
from typing import List, Optional

from sqlalchemy import (
    JSON,
    BigInteger,
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY, UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base


# Helper function to generate UUID defaults
def generate_uuid() -> uuid.UUID:
    """
    Generate a new UUID for use as default.

    Returns:
        A new UUID4 instance
    """
    return uuid.uuid4()


class GateType(str, PyEnum):
    """Types of flow cytometry gates."""
    POLYGON = "polygon"
    RECTANGLE = "rectangle"
    QUADRANT = "quadrant"


class BooleanOperator(str, PyEnum):
    """Boolean operators for combining gates."""
    AND = "AND"
    OR = "OR"
    NOT = "NOT"


class FlowCytometryDataset(Base):
    """Flow cytometry dataset metadata linked to a Dataset."""

    __tablename__ = "flow_cytometry_datasets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    dataset_id = Column(UUID(as_uuid=True), ForeignKey("datasets.id"), nullable=False, index=True)

    # File storage - events stored in Parquet format for performance
    events_parquet_path = Column(String(500), nullable=False)
    file_size_bytes = Column(BigInteger, nullable=True)
    
    # Dataset statistics
    n_events = Column(Integer, nullable=True)
    n_parameters = Column(Integer, nullable=True)
    
    # Acquisition metadata
    acquisition_date = Column(DateTime(timezone=True), nullable=True)
    cytometer_model = Column(String(200), nullable=True)
    cytometer_serial = Column(String(100), nullable=True)
    acquisition_software = Column(String(100), nullable=True)
    acquisition_settings = Column(JSON, nullable=True)  # Voltage settings, compensation matrix, etc.
    
    # Sample information
    sample_id = Column(String(200), nullable=True)
    sample_volume_ul = Column(Float, nullable=True)
    dilution_factor = Column(Float, nullable=True)
    staining_protocol = Column(Text, nullable=True)
    
    # Processing status
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    
    # Timestamps
    ingested_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    processed_at = Column(DateTime(timezone=True), nullable=True)

    # Relationships
    dataset: Mapped["Dataset"] = relationship("Dataset", foreign_keys=[dataset_id])
    parameters: Mapped[List["FlowCytometryParameter"]] = relationship(
        "FlowCytometryParameter", back_populates="flow_dataset", cascade="all, delete-orphan"
    )
    gates: Mapped[List["FlowCytometryGate"]] = relationship(
        "FlowCytometryGate", back_populates="flow_dataset", cascade="all, delete-orphan"
    )
    populations: Mapped[List["FlowCytometryPopulation"]] = relationship(
        "FlowCytometryPopulation", back_populates="flow_dataset", cascade="all, delete-orphan"
    )


class FlowCytometryParameter(Base):
    """Channel/parameter metadata for a flow cytometry dataset."""

    __tablename__ = "flow_cytometry_parameters"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("flow_cytometry_datasets.id"), nullable=False, index=True
    )
    
    parameter_index = Column(Integer, nullable=False)  # 0-based index in the data
    parameter_name = Column(String(100), nullable=False)  # e.g., "FSC-A", "CD3-FITC"
    parameter_short_name = Column(String(50), nullable=True)  # e.g., "FSC", "CD3"
    
    # Channel configuration
    detector = Column(String(100), nullable=True)  # PMT, APD, etc.
    excitation_wavelength = Column(Integer, nullable=True)  # nm
    emission_wavelength = Column(Integer, nullable=True)  # nm
    fluorophore = Column(String(100), nullable=True)  # FITC, PE, APC, etc.
    
    # Data range and scaling
    data_range = Column(Integer, nullable=True)  # e.g., 262144 for 18-bit
    amplifier_gain = Column(Float, nullable=True)
    voltage = Column(Float, nullable=True)
    
    # Statistics
    min_value = Column(Float, nullable=True)
    max_value = Column(Float, nullable=True)
    
    # Relationships
    flow_dataset: Mapped["FlowCytometryDataset"] = relationship(
        "FlowCytometryDataset", foreign_keys=[flow_dataset_id], back_populates="parameters"
    )

    __table_args__ = (
        UniqueConstraint("flow_dataset_id", "parameter_index", name="uq_flow_parameter_index"),
        Index("ix_flow_parameters_dataset_name", "flow_dataset_id", "parameter_name"),
    )


class FlowCytometryGate(Base):
    """Gate definition for flow cytometry analysis."""

    __tablename__ = "flow_cytometry_gates"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("flow_cytometry_datasets.id"), nullable=False, index=True
    )
    
    gate_name = Column(String(200), nullable=False)
    gate_type = Column(String(20), nullable=False)  # "polygon", "rectangle", "quadrant"
    
    # Gate geometry - stored as JSON for flexibility
    # polygon: {"vertices": [[x1, y1], [x2, y2], ...]}
    # rectangle: {"x_min": x1, "x_max": x2, "y_min": y1, "y_max": y2}
    # quadrant: {"x_threshold": x, "y_threshold": y, "quadrant": "LL|LR|UL|UR"}
    gate_definition = Column(JSON, nullable=False)
    
    # Parameters this gate operates on
    x_parameter_id = Column(UUID(as_uuid=True), ForeignKey("flow_cytometry_parameters.id"), nullable=False)
    y_parameter_id = Column(UUID(as_uuid=True), ForeignKey("flow_cytometry_parameters.id"), nullable=True)  # Null for 1D gates
    
    # Boolean logic for combining gates
    boolean_operator = Column(String(10), nullable=True)  # "AND", "OR", "NOT"
    operand_gate_ids = Column(ARRAY(UUID(as_uuid=True)), nullable=True)  # Gates to combine with this one
    
    # Gate metadata
    parent_gate_id = Column(UUID(as_uuid=True), ForeignKey("flow_cytometry_gates.id"), nullable=True)
    is_active = Column(Boolean, nullable=False, default=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationships
    flow_dataset: Mapped["FlowCytometryDataset"] = relationship(
        "FlowCytometryDataset", foreign_keys=[flow_dataset_id], back_populates="gates"
    )
    x_parameter: Mapped["FlowCytometryParameter"] = relationship(
        "FlowCytometryParameter", foreign_keys=[x_parameter_id]
    )
    y_parameter: Mapped[Optional["FlowCytometryParameter"]] = relationship(
        "FlowCytometryParameter", foreign_keys=[y_parameter_id]
    )
    parent_gate: Mapped[Optional["FlowCytometryGate"]] = relationship(
        "FlowCytometryGate", foreign_keys=[parent_gate_id], remote_side=[id]
    )
    child_gates: Mapped[List["FlowCytometryGate"]] = relationship(
        "FlowCytometryGate", foreign_keys=[parent_gate_id], cascade="all, delete-orphan"
    )
    populations: Mapped[List["FlowCytometryPopulation"]] = relationship(
        "FlowCytometryPopulation", back_populates="gate", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_flow_gates_dataset_name", "flow_dataset_id", "gate_name"),
        Index("ix_flow_gates_parent", "parent_gate_id"),
    )


class FlowCytometryPopulation(Base):
    """Statistics for a gated population in flow cytometry analysis."""

    __tablename__ = "flow_cytometry_populations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    flow_dataset_id = Column(
        UUID(as_uuid=True), ForeignKey("flow_cytometry_datasets.id"), nullable=False, index=True
    )
    gate_id = Column(UUID(as_uuid=True), ForeignKey("flow_cytometry_gates.id"), nullable=False, index=True)
    
    # Population statistics
    event_count = Column(Integer, nullable=False)
    parent_event_count = Column(Integer, nullable=True)  # Events in parent gate
    percentage_of_parent = Column(Float, nullable=True)
    percentage_of_total = Column(Float, nullable=True)
    
    # Population characteristics
    population_name = Column(String(200), nullable=True)
    cell_type = Column(String(200), nullable=True)
    phenotype = Column(String(500), nullable=True)  # e.g., "CD3+ CD4+ CD8-"
    
    # Statistical measures for each parameter
    parameter_statistics = Column(JSON, nullable=True)  # {"FSC-A": {"mean": 1000, "median": 950, "std": 200}, ...}
    
    # Analysis metadata
    analysis_software = Column(String(100), nullable=True)
    analysis_version = Column(String(50), nullable=True)
    analyzed_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationships
    flow_dataset: Mapped["FlowCytometryDataset"] = relationship(
        "FlowCytometryDataset", foreign_keys=[flow_dataset_id], back_populates="populations"
    )
    gate: Mapped["FlowCytometryGate"] = relationship(
        "FlowCytometryGate", foreign_keys=[gate_id], back_populates="populations"
    )

    __table_args__ = (
        UniqueConstraint("flow_dataset_id", "gate_id", name="uq_flow_population_gate"),
        Index("ix_flow_populations_dataset_gate", "flow_dataset_id", "gate_id"),
    )


# Export all models
__all__ = [
    "FlowCytometryDataset",
    "FlowCytometryParameter", 
    "FlowCytometryGate",
    "FlowCytometryPopulation",
    "GateType",
    "BooleanOperator",
]
