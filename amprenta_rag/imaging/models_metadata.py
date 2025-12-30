"""Imaging metadata database models for microscope instruments and acquisition settings."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, JSON, String, Text, Index
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.imaging.models import MicroscopyImage
    from amprenta_rag.models.chemistry import HTSPlate


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class Microscope(Base):
    """Microscope instrument registry."""

    __tablename__ = "microscopes"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(200), nullable=False)  # "Nikon Ti2-E", "Zeiss Axio Observer"
    manufacturer = Column(String(100), nullable=False)  # "Nikon", "Zeiss", "Leica"
    model = Column(String(100), nullable=False)
    serial_number = Column(String(100), unique=True, nullable=True)
    facility_location = Column(String(200), nullable=True)  # "Building A, Room 101"
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    objectives: Mapped[List["Objective"]] = relationship(
        back_populates="microscope", cascade="all, delete-orphan"
    )
    light_sources: Mapped[List["LightSource"]] = relationship(
        back_populates="microscope", cascade="all, delete-orphan"
    )
    channel_configs: Mapped[List["ChannelConfig"]] = relationship(
        back_populates="microscope", cascade="all, delete-orphan"
    )

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_microscopes_manufacturer_model", "manufacturer", "model"),
        Index("idx_microscopes_active", "is_active"),
    )


class Objective(Base):
    """Objective lens specifications."""

    __tablename__ = "objectives"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    microscope_id = Column(UUID(as_uuid=True), ForeignKey("microscopes.id"), nullable=True, index=True)  # Can be shared
    name = Column(String(100), nullable=False)  # "Plan Apo 20x"
    magnification = Column(Float, nullable=False)  # 20.0
    numerical_aperture = Column(Float, nullable=False)  # 0.75
    immersion = Column(String(50), nullable=False)  # "air", "oil", "water", "glycerol"
    working_distance_mm = Column(Float, nullable=True)
    correction = Column(String(50), nullable=True)  # "Plan", "Fluor", "Apo"
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    microscope: Mapped[Optional["Microscope"]] = relationship(back_populates="objectives")

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_objectives_magnification", "magnification"),
        Index("idx_objectives_active", "is_active"),
    )


class LightSource(Base):
    """Laser/LED light source."""

    __tablename__ = "light_sources"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    microscope_id = Column(UUID(as_uuid=True), ForeignKey("microscopes.id"), nullable=False, index=True)
    name = Column(String(100), nullable=False)  # "488nm Argon Laser"
    source_type = Column(String(50), nullable=False)  # "laser", "led", "halogen", "arc"
    wavelength_nm = Column(Integer, nullable=True)  # 488
    max_power_mw = Column(Float, nullable=True)
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    microscope: Mapped["Microscope"] = relationship(back_populates="light_sources")
    channel_configs: Mapped[List["ChannelConfig"]] = relationship(back_populates="light_source")

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_light_sources_wavelength", "wavelength_nm"),
        Index("idx_light_sources_type", "source_type"),
    )


class FilterSet(Base):
    """Excitation/emission filter configuration."""

    __tablename__ = "filter_sets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    name = Column(String(100), nullable=False)  # "FITC", "DAPI", "Texas Red"
    excitation_center_nm = Column(Integer, nullable=False)
    excitation_bandwidth_nm = Column(Integer, nullable=False)
    emission_center_nm = Column(Integer, nullable=False)
    emission_bandwidth_nm = Column(Integer, nullable=False)
    dichroic_cutoff_nm = Column(Integer, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    channel_configs: Mapped[List["ChannelConfig"]] = relationship(back_populates="filter_set")

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_filter_sets_excitation", "excitation_center_nm"),
        Index("idx_filter_sets_emission", "emission_center_nm"),
    )


class ChannelConfig(Base):
    """Channel-fluorophore-filter mapping."""

    __tablename__ = "channel_configs"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    microscope_id = Column(UUID(as_uuid=True), ForeignKey("microscopes.id"), nullable=False, index=True)
    channel_name = Column(String(50), nullable=False)  # "DAPI", "GFP", "RFP"
    fluorophore = Column(String(100), nullable=True)  # "Hoechst 33342"
    light_source_id = Column(UUID(as_uuid=True), ForeignKey("light_sources.id"), nullable=True, index=True)
    filter_set_id = Column(UUID(as_uuid=True), ForeignKey("filter_sets.id"), nullable=True, index=True)
    default_exposure_ms = Column(Float, nullable=True)
    default_gain = Column(Float, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    microscope: Mapped["Microscope"] = relationship(back_populates="channel_configs")
    light_source: Mapped[Optional["LightSource"]] = relationship(back_populates="channel_configs")
    filter_set: Mapped[Optional["FilterSet"]] = relationship(back_populates="channel_configs")

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_channel_configs_microscope_channel", "microscope_id", "channel_name"),
    )


class AcquisitionSettings(Base):
    """Per-image acquisition parameters."""

    __tablename__ = "acquisition_settings"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    # FK back to MicroscopyImage is on MicroscopyImage side (1:1)
    exposure_ms = Column(Float, nullable=False)
    gain = Column(Float, nullable=True)
    laser_power_percent = Column(Float, nullable=True)
    binning = Column(Integer, default=1, nullable=False)  # 1x1, 2x2, etc.
    z_position_um = Column(Float, nullable=True)
    autofocus_offset_um = Column(Float, nullable=True)
    temperature_celsius = Column(Float, nullable=True)
    raw_settings = Column(JSON, nullable=True)  # Vendor-specific settings
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships (1:1 with MicroscopyImage)
    microscopy_image: Mapped[Optional["MicroscopyImage"]] = relationship(back_populates="acquisition_settings")


class ImageFileSet(Base):
    """Batch import tracking for multi-file imaging datasets."""

    __tablename__ = "image_file_sets"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    plate_id = Column(UUID(as_uuid=True), ForeignKey("hts_plates.id"), nullable=True, index=True)
    vendor = Column(String(100), nullable=False)  # "opera", "imagexpress", "cellvoyager"
    import_path = Column(String(500), nullable=False)
    file_count = Column(Integer, nullable=False)
    image_count = Column(Integer, default=0, nullable=False)
    import_status = Column(String(50), default="pending", nullable=False)  # pending/importing/completed/failed
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    plate: Mapped[Optional["HTSPlate"]] = relationship()

    # Indexes for efficient queries
    __table_args__ = (
        Index("idx_image_file_sets_status", "import_status"),
        Index("idx_image_file_sets_vendor", "vendor"),
    )
