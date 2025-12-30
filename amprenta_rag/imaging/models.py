"""Imaging and microscopy database models."""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional, TYPE_CHECKING

from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, JSON, String, Text, Index
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, relationship

from amprenta_rag.database.base import Base

if TYPE_CHECKING:
    from amprenta_rag.models.chemistry import HTSWell
    from amprenta_rag.imaging.models_metadata import Objective, AcquisitionSettings


def generate_uuid() -> uuid.UUID:
    return uuid.uuid4()


class MicroscopyImage(Base):
    """Microscopy image from HTS plate imaging."""

    __tablename__ = "microscopy_images"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    well_id = Column(UUID(as_uuid=True), ForeignKey("hts_wells.id"), nullable=True, index=True)  # P0 fix: nullable for standalone OME-TIFF
    channel = Column(String(50), nullable=False)  # e.g., "DAPI", "GFP", "RFP", "Brightfield"
    z_slice = Column(Integer, nullable=True, default=0)  # Z-stack slice index
    timepoint = Column(Integer, nullable=True, default=0)  # Time series index
    
    # Metadata relationships (new)
    objective_id = Column(UUID(as_uuid=True), ForeignKey("objectives.id"), nullable=True, index=True)
    acquisition_settings_id = Column(UUID(as_uuid=True), ForeignKey("acquisition_settings.id"), nullable=True, index=True)
    
    # Image dimensions and metadata
    width = Column(Integer, nullable=False)
    height = Column(Integer, nullable=False)
    bit_depth = Column(Integer, nullable=False, default=16)
    pixel_size_um = Column(Float, nullable=True)  # Microns per pixel
    
    # File storage paths
    image_path = Column(String(500), nullable=False)  # Path to image file
    thumbnail_path = Column(String(500), nullable=True)  # Path to thumbnail
    
    # Image metadata
    image_metadata = Column(JSON, nullable=True)  # Acquisition settings, exposure, etc.
    ome_metadata = Column(JSON, nullable=True)  # P2 fix: OME metadata for queryability
    
    # Quality metrics
    focus_score = Column(Float, nullable=True)  # Image focus quality
    signal_noise_ratio = Column(Float, nullable=True)
    
    # Timestamps
    acquired_at = Column(DateTime, nullable=True)  # When image was acquired
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Relationships
    well: Mapped[Optional["HTSWell"]] = relationship()  # Optional for standalone OME-TIFF
    objective: Mapped[Optional["Objective"]] = relationship("amprenta_rag.imaging.models_metadata.Objective")
    acquisition_settings: Mapped[Optional["AcquisitionSettings"]] = relationship("amprenta_rag.imaging.models_metadata.AcquisitionSettings", back_populates="microscopy_image")
    segmentations: Mapped[List["CellSegmentation"]] = relationship(
        back_populates="image", cascade="all, delete-orphan"
    )

    # Index for efficient plate-level queries
    __table_args__ = (
        Index("idx_microscopy_images_well_channel", "well_id", "channel"),
        Index("idx_microscopy_images_well_z_time", "well_id", "z_slice", "timepoint"),
    )


class CellSegmentation(Base):
    """Cell segmentation results from image analysis."""

    __tablename__ = "cell_segmentations"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    image_id = Column(UUID(as_uuid=True), ForeignKey("microscopy_images.id"), nullable=False, index=True)
    model_name = Column(String(200), nullable=False)  # e.g., "cellpose", "stardist", "custom_unet"
    model_version = Column(String(100), nullable=True)
    
    # Segmentation results
    cell_count = Column(Integer, nullable=False, default=0)
    mask_path = Column(String(500), nullable=False)  # Path to segmentation mask (NPY format)
    
    # Segmentation parameters
    parameters = Column(JSON, nullable=True)  # Model parameters used
    
    # Quality metrics
    confidence_score = Column(Float, nullable=True)  # Overall segmentation confidence
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)

    # Relationships
    image: Mapped["MicroscopyImage"] = relationship(back_populates="segmentations")
    features: Mapped[List["CellFeature"]] = relationship(
        back_populates="segmentation", cascade="all, delete-orphan"
    )


class CellFeature(Base):
    """Extracted features from segmented cells."""

    __tablename__ = "cell_features"

    id = Column(UUID(as_uuid=True), primary_key=True, default=generate_uuid)
    segmentation_id = Column(UUID(as_uuid=True), ForeignKey("cell_segmentations.id"), nullable=False, index=True)
    cell_id = Column(Integer, nullable=False)  # Cell ID within the segmentation mask
    
    # Morphological features
    area = Column(Float, nullable=True)  # Cell area in pixels
    perimeter = Column(Float, nullable=True)  # Cell perimeter in pixels
    circularity = Column(Float, nullable=True)  # Shape circularity (0-1)
    eccentricity = Column(Float, nullable=True)  # Shape eccentricity (0-1)
    solidity = Column(Float, nullable=True)  # Convex hull solidity
    
    # Position features
    centroid_x = Column(Float, nullable=True)  # X coordinate of centroid
    centroid_y = Column(Float, nullable=True)  # Y coordinate of centroid
    
    # Intensity features (per channel)
    intensity_features = Column(JSON, nullable=True)  # Mean, median, std, min, max per channel
    
    # Texture features
    texture_features = Column(JSON, nullable=True)  # GLCM, LBP, etc.
    
    # Additional custom features
    custom_features = Column(JSON, nullable=True)  # Experiment-specific features
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    segmentation: Mapped["CellSegmentation"] = relationship(back_populates="features")

    # Index for efficient queries
    __table_args__ = (
        Index("idx_cell_features_segmentation_cell", "segmentation_id", "cell_id"),
    )
