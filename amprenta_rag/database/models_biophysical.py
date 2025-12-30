"""
Biophysical Assay database models for the multi-omics platform.

These models store SPR (Surface Plasmon Resonance), MST (Microscale Thermophoresis),
and DSC (Differential Scanning Calorimetry) experimental data and results.
Time-series data is stored in ARRAY columns for performance.
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


class SPRExperimentType(str, PyEnum):
    """Types of SPR experiments."""
    KINETICS = "kinetics"
    AFFINITY = "affinity"
    SCREENING = "screening"


class SPRExperiment(Base):
    """SPR (Surface Plasmon Resonance) experiment metadata and kinetic parameters."""

    __tablename__ = "spr_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    
    # Experiment identification
    experiment_name = Column(String(200), nullable=False, index=True)
    experiment_type = Column(String(20), nullable=False, default=SPRExperimentType.KINETICS)
    
    # Target and compound relationships
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)  # Fallback if no target_id
    
    # Experimental conditions
    analyte_name = Column(String(200), nullable=True)  # What was injected
    ligand_name = Column(String(200), nullable=True)   # What was immobilized
    analyte_concentration = Column(Float, nullable=True)
    ligand_density = Column(Float, nullable=True)  # RU immobilized
    buffer_composition = Column(String(500), nullable=True)
    temperature_celsius = Column(Float, nullable=True, default=25.0)
    flow_rate = Column(Float, nullable=True)  # μL/min
    
    # Units
    response_units = Column(String(50), nullable=False, default="RU")
    concentration_units = Column(String(20), nullable=False, default="nM")
    
    # Kinetic parameters (calculated from sensorgrams)
    ka = Column(Float, nullable=True)  # Association rate constant (1/Ms)
    kd = Column(Float, nullable=True)  # Dissociation rate constant (1/s)
    kd_equilibrium = Column(Float, nullable=True)  # Equilibrium dissociation constant (M)
    rmax = Column(Float, nullable=True)  # Maximum binding response (RU)
    chi_squared = Column(Float, nullable=True)  # Goodness of fit
    
    # Data quality metrics
    baseline_drift = Column(Float, nullable=True)  # RU/s
    bulk_refractive_index = Column(Float, nullable=True)
    reference_subtracted = Column(Boolean, nullable=False, default=True)
    
    # Experiment metadata
    instrument_model = Column(String(100), nullable=True)  # e.g., "Biacore T200"
    chip_type = Column(String(50), nullable=True)  # e.g., "CM5", "NTA"
    cycle_number = Column(Integer, nullable=True)
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    
    # Processing status
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    
    # Timestamps
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    
    # Relationships
    compound: Mapped[Optional["Compound"]] = relationship("Compound", foreign_keys=[compound_id])
    target: Mapped[Optional["ProteinStructure"]] = relationship("ProteinStructure", foreign_keys=[target_id])
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    sensorgrams: Mapped[List["SPRSensorgram"]] = relationship(
        "SPRSensorgram", back_populates="experiment", cascade="all, delete-orphan"
    )
    
    __table_args__ = (
        Index("ix_spr_experiments_compound_target", "compound_id", "target_id"),
        Index("ix_spr_experiments_replicate", "replicate_group_id", "replicate_number"),
    )


class SPRSensorgram(Base):
    """SPR sensorgram time-series data (max 50,000 data points for performance)."""

    __tablename__ = "spr_sensorgrams"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("spr_experiments.id"), nullable=False, index=True)
    
    # Sensorgram identification
    cycle_phase = Column(String(50), nullable=False)  # "association", "dissociation", "regeneration"
    flow_cell = Column(Integer, nullable=True)  # Flow cell number (1-4)
    
    # Time-series data (max 50,000 points for performance)
    time_seconds = Column(ARRAY(Float), nullable=True)  # Max 50,000 points
    response_values = Column(ARRAY(Float), nullable=True)  # Response in RU
    
    # Data processing flags
    baseline_corrected = Column(Boolean, nullable=False, default=False)
    reference_subtracted = Column(Boolean, nullable=False, default=False)
    bulk_shift_corrected = Column(Boolean, nullable=False, default=False)
    
    # Data quality
    data_points_count = Column(Integer, nullable=True)
    noise_level = Column(Float, nullable=True)  # RMS noise in RU
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationships
    experiment: Mapped["SPRExperiment"] = relationship("SPRExperiment", back_populates="sensorgrams")
    
    __table_args__ = (
        Index("ix_spr_sensorgrams_experiment_phase", "experiment_id", "cycle_phase"),
    )


class MSTExperiment(Base):
    """MST (Microscale Thermophoresis) experiment metadata and affinity measurements."""

    __tablename__ = "mst_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    
    # Experiment identification
    experiment_name = Column(String(200), nullable=False, index=True)
    
    # Target and compound relationships
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)  # Fallback if no target_id
    
    # Experimental conditions
    target_concentration = Column(Float, nullable=True)  # Constant target concentration
    ligand_name = Column(String(200), nullable=True)    # Titrated ligand
    buffer_composition = Column(String(500), nullable=True)
    temperature_celsius = Column(Float, nullable=True, default=25.0)
    excitation_power = Column(Float, nullable=True)  # % LED power
    mst_power = Column(Float, nullable=True)  # % MST power
    
    # Units
    response_units = Column(String(50), nullable=False, default="‰")  # per mille
    concentration_units = Column(String(20), nullable=False, default="nM")
    
    # Binding affinity results
    kd_value = Column(Float, nullable=True)  # Dissociation constant
    kd_error = Column(Float, nullable=True)  # Standard error of KD
    hill_coefficient = Column(Float, nullable=True)
    binding_amplitude = Column(Float, nullable=True)  # Change in thermophoresis
    r_squared = Column(Float, nullable=True)  # Goodness of fit
    
    # Data quality metrics
    signal_to_noise = Column(Float, nullable=True)
    aggregation_detected = Column(Boolean, nullable=False, default=False)
    bleaching_rate = Column(Float, nullable=True)  # % per second
    
    # Experiment metadata
    instrument_model = Column(String(100), nullable=True)  # e.g., "NanoTemper Monolith NT.115"
    capillary_type = Column(String(50), nullable=True)  # Standard, Premium, etc.
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    
    # Processing status
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    
    # Timestamps
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    
    # Relationships
    compound: Mapped[Optional["Compound"]] = relationship("Compound", foreign_keys=[compound_id])
    target: Mapped[Optional["ProteinStructure"]] = relationship("ProteinStructure", foreign_keys=[target_id])
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    dose_responses: Mapped[List["MSTDoseResponse"]] = relationship(
        "MSTDoseResponse", back_populates="experiment", cascade="all, delete-orphan"
    )
    
    __table_args__ = (
        Index("ix_mst_experiments_compound_target", "compound_id", "target_id"),
        Index("ix_mst_experiments_replicate", "replicate_group_id", "replicate_number"),
    )


class MSTDoseResponse(Base):
    """MST dose-response curve data points."""

    __tablename__ = "mst_dose_responses"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("mst_experiments.id"), nullable=False, index=True)
    
    # Dose-response data
    ligand_concentration = Column(Float, nullable=False)  # Concentration of titrated ligand
    thermophoresis_response = Column(Float, nullable=False)  # Normalized thermophoresis signal
    thermophoresis_error = Column(Float, nullable=True)  # Standard error
    
    # Individual measurement data
    initial_fluorescence = Column(Float, nullable=True)  # F_initial
    thermophoresis_signal = Column(Float, nullable=True)  # Raw thermophoresis
    
    # Data quality for this point
    capillary_number = Column(Integer, nullable=True)
    measurement_valid = Column(Boolean, nullable=False, default=True)
    outlier_flag = Column(Boolean, nullable=False, default=False)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationships
    experiment: Mapped["MSTExperiment"] = relationship("MSTExperiment", back_populates="dose_responses")
    
    __table_args__ = (
        Index("ix_mst_dose_responses_experiment_conc", "experiment_id", "ligand_concentration"),
    )


class DSCExperiment(Base):
    """DSC (Differential Scanning Calorimetry) experiment metadata and thermal parameters."""

    __tablename__ = "dsc_experiments"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    
    # Experiment identification
    experiment_name = Column(String(200), nullable=False, index=True)
    
    # Target and compound relationships
    compound_id = Column(UUID(as_uuid=True), ForeignKey("compounds.id"), nullable=True, index=True)
    target_id = Column(UUID(as_uuid=True), ForeignKey("protein_structures.id"), nullable=True, index=True)
    target_name = Column(String(200), nullable=True)  # Fallback if no target_id
    
    # Experimental conditions
    protein_concentration = Column(Float, nullable=True)  # mg/mL or μM
    ligand_concentration = Column(Float, nullable=True)  # For binding studies
    buffer_composition = Column(String(500), nullable=True)
    scan_rate = Column(Float, nullable=True, default=1.0)  # °C/min
    start_temperature = Column(Float, nullable=True, default=20.0)  # °C
    end_temperature = Column(Float, nullable=True, default=100.0)  # °C
    
    # Units
    response_units = Column(String(50), nullable=False, default="kcal/mol/°C")
    concentration_units = Column(String(20), nullable=False, default="μM")
    
    # Thermal parameters (calculated from thermogram)
    tm_value = Column(Float, nullable=True)  # Melting temperature (°C)
    tm_error = Column(Float, nullable=True)  # Standard error of Tm
    delta_h = Column(Float, nullable=True)  # Enthalpy change (kcal/mol)
    delta_h_error = Column(Float, nullable=True)
    delta_cp = Column(Float, nullable=True)  # Heat capacity change
    cooperativity = Column(Float, nullable=True)  # Cooperative unit size
    
    # Ligand binding parameters (if applicable)
    kd_thermal = Column(Float, nullable=True)  # KD from thermal shift
    delta_tm = Column(Float, nullable=True)  # Thermal shift (°C)
    
    # Data quality metrics
    baseline_quality = Column(Float, nullable=True)  # R² of baseline fit
    peak_symmetry = Column(Float, nullable=True)
    signal_to_noise = Column(Float, nullable=True)
    
    # Experiment metadata
    instrument_model = Column(String(100), nullable=True)  # e.g., "MicroCal PEAQ-DSC"
    cell_volume = Column(Float, nullable=True)  # μL
    reference_buffer = Column(String(500), nullable=True)
    replicate_number = Column(Integer, nullable=False, default=1)
    replicate_group_id = Column(UUID(as_uuid=True), nullable=True, index=True)
    
    # Processing status
    processing_status = Column(String(50), nullable=False, default="pending")
    processing_log = Column(Text, nullable=True)
    
    # Timestamps
    experiment_date = Column(DateTime(timezone=True), nullable=True)
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    updated_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), onupdate=lambda: datetime.now(timezone.utc), nullable=False)
    created_by_id = Column(UUID(as_uuid=True), ForeignKey("users.id"), nullable=True)
    
    # Relationships
    compound: Mapped[Optional["Compound"]] = relationship("Compound", foreign_keys=[compound_id])
    target: Mapped[Optional["ProteinStructure"]] = relationship("ProteinStructure", foreign_keys=[target_id])
    created_by: Mapped[Optional["User"]] = relationship("User", foreign_keys=[created_by_id])
    scans: Mapped[List["DSCScan"]] = relationship(
        "DSCScan", back_populates="experiment", cascade="all, delete-orphan"
    )
    
    __table_args__ = (
        Index("ix_dsc_experiments_compound_target", "compound_id", "target_id"),
        Index("ix_dsc_experiments_replicate", "replicate_group_id", "replicate_number"),
    )


class DSCScan(Base):
    """DSC thermogram scan data (max 50,000 data points for performance)."""

    __tablename__ = "dsc_scans"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey("dsc_experiments.id"), nullable=False, index=True)
    
    # Scan identification
    scan_type = Column(String(50), nullable=False)  # "sample", "reference", "buffer"
    scan_number = Column(Integer, nullable=True)
    
    # Thermogram data (max 50,000 points for performance)
    temperature_celsius = Column(ARRAY(Float), nullable=True)  # Max 50,000 points
    heat_capacity = Column(ARRAY(Float), nullable=True)  # kcal/mol/°C
    
    # Data processing flags
    baseline_subtracted = Column(Boolean, nullable=False, default=False)
    reference_subtracted = Column(Boolean, nullable=False, default=False)
    normalized = Column(Boolean, nullable=False, default=False)
    
    # Data quality
    data_points_count = Column(Integer, nullable=True)
    noise_level = Column(Float, nullable=True)  # RMS noise
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationships
    experiment: Mapped["DSCExperiment"] = relationship("DSCExperiment", back_populates="scans")
    
    __table_args__ = (
        Index("ix_dsc_scans_experiment_type", "experiment_id", "scan_type"),
    )


# Export all models
__all__ = [
    "SPRExperimentType",
    "SPRExperiment",
    "SPRSensorgram", 
    "MSTExperiment",
    "MSTDoseResponse",
    "DSCExperiment",
    "DSCScan",
]
