"""Pydantic schemas for imaging API endpoints."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

from pydantic import BaseModel, Field

# Define BaseSchema locally to avoid circular imports
from pydantic import BaseModel, ConfigDict

class BaseSchema(BaseModel):
    """Base schema with common configuration."""

    model_config = ConfigDict(
        from_attributes=True,
        json_encoders={
            UUID: str,
            datetime: lambda v: v.isoformat(),
        }
    )


# ============================================================================
# OME-TIFF Import schemas
# ============================================================================


class OMETiffImportRequest(BaseSchema):
    """Request schema for OME-TIFF import."""
    
    well_id: Optional[UUID] = Field(None, description="Well ID to associate with images")
    plate_id: Optional[UUID] = Field(None, description="Plate ID for standalone import")
    channel_mapping: Optional[Dict[str, str]] = Field(
        None, 
        description="Channel index to name mapping (e.g., {'0': 'DAPI', '1': 'GFP'})"
    )


class OMETiffImportResponse(BaseSchema):
    """Response schema for OME-TIFF import."""
    
    image_id: UUID = Field(description="Primary image ID")
    filename: str = Field(description="Original filename")
    dimensions: Dict[str, int] = Field(description="Image dimensions (x, y, z, c, t)")
    channels: List[str] = Field(description="Channel names")
    instrument: Optional[str] = Field(None, description="Microscope instrument name")
    pixel_size_um: Optional[float] = Field(None, description="Pixel size in microns")
    ome_metadata: Dict[str, Any] = Field(description="Parsed OME metadata")


# ============================================================================
# Batch Import schemas
# ============================================================================


class BatchImportRequest(BaseSchema):
    """Request schema for batch vendor import."""
    
    vendor: Optional[str] = Field(None, description="Vendor format (auto-detect if None)")
    plate_id: Optional[UUID] = Field(None, description="Existing plate ID")
    campaign_id: Optional[UUID] = Field(None, description="Campaign to associate with")
    create_missing_wells: bool = Field(True, description="Create wells if they don't exist")


class BatchImportResponse(BaseSchema):
    """Response schema for batch import."""
    
    fileset_id: UUID = Field(description="Import job identifier")
    status: str = Field(description="Import status (pending, importing, completed, failed)")
    vendor: str = Field(description="Detected or specified vendor")
    total_images: int = Field(description="Total images to import")
    imported_count: int = Field(description="Images imported so far")
    errors: List[str] = Field(description="Import error messages")


class ImportStatusResponse(BaseSchema):
    """Response schema for import status."""
    
    fileset_id: UUID = Field(description="Import job identifier")
    status: str = Field(description="Current import status")
    progress_percent: float = Field(description="Import progress (0-100)")
    total_images: int = Field(description="Total images to import")
    imported_count: int = Field(description="Images imported successfully")
    failed_count: int = Field(description="Images that failed import")
    errors: List[str] = Field(description="Error messages")
    warnings: List[str] = Field(description="Warning messages")
    estimated_completion: Optional[datetime] = Field(None, description="Estimated completion time")


# ============================================================================
# Instrument Management schemas
# ============================================================================


class MicroscopeCreate(BaseSchema):
    """Request schema for creating a microscope."""
    
    name: str = Field(..., min_length=1, max_length=200, description="Microscope name")
    manufacturer: str = Field(..., min_length=1, max_length=100, description="Manufacturer")
    model: str = Field(..., min_length=1, max_length=100, description="Model")
    serial_number: Optional[str] = Field(None, max_length=100, description="Serial number")
    facility_location: Optional[str] = Field(None, max_length=200, description="Location")


class MicroscopeResponse(BaseSchema):
    """Response schema for microscope."""
    
    id: UUID = Field(description="Microscope ID")
    name: str = Field(description="Microscope name")
    manufacturer: str = Field(description="Manufacturer")
    model: str = Field(description="Model")
    serial_number: Optional[str] = Field(description="Serial number")
    facility_location: Optional[str] = Field(description="Facility location")
    is_active: bool = Field(description="Whether microscope is active")
    objectives: List["ObjectiveResponse"] = Field(description="Associated objectives")
    channels: List["ChannelConfigResponse"] = Field(description="Channel configurations")
    created_at: datetime = Field(description="Creation timestamp")


class ObjectiveCreate(BaseSchema):
    """Request schema for creating an objective."""
    
    microscope_id: Optional[UUID] = Field(None, description="Microscope ID (can be shared)")
    name: str = Field(..., min_length=1, max_length=100, description="Objective name")
    magnification: float = Field(..., gt=0, description="Magnification (e.g., 20.0)")
    numerical_aperture: float = Field(..., gt=0, le=2.0, description="Numerical aperture")
    immersion: str = Field(..., description="Immersion medium (air, oil, water, glycerol)")
    working_distance_mm: Optional[float] = Field(None, gt=0, description="Working distance in mm")
    correction: Optional[str] = Field(None, description="Correction type (Plan, Fluor, Apo)")


class ObjectiveResponse(BaseSchema):
    """Response schema for objective."""
    
    id: UUID = Field(description="Objective ID")
    microscope_id: Optional[UUID] = Field(description="Associated microscope ID")
    name: str = Field(description="Objective name")
    magnification: float = Field(description="Magnification")
    numerical_aperture: float = Field(description="Numerical aperture")
    immersion: str = Field(description="Immersion medium")
    working_distance_mm: Optional[float] = Field(description="Working distance in mm")
    correction: Optional[str] = Field(description="Correction type")
    is_active: bool = Field(description="Whether objective is active")


class ChannelConfigCreate(BaseSchema):
    """Request schema for creating a channel configuration."""
    
    microscope_id: UUID = Field(..., description="Microscope ID")
    channel_name: str = Field(..., min_length=1, max_length=50, description="Channel name")
    fluorophore: Optional[str] = Field(None, max_length=100, description="Fluorophore name")
    light_source_id: Optional[UUID] = Field(None, description="Light source ID")
    filter_set_id: Optional[UUID] = Field(None, description="Filter set ID")
    default_exposure_ms: Optional[float] = Field(None, gt=0, description="Default exposure time")
    default_gain: Optional[float] = Field(None, ge=0, description="Default gain")


class ChannelConfigResponse(BaseSchema):
    """Response schema for channel configuration."""
    
    id: UUID = Field(description="Channel config ID")
    microscope_id: UUID = Field(description="Microscope ID")
    channel_name: str = Field(description="Channel name")
    fluorophore: Optional[str] = Field(description="Fluorophore name")
    excitation_wavelength_nm: Optional[int] = Field(description="Excitation wavelength")
    emission_wavelength_nm: Optional[int] = Field(description="Emission wavelength")
    default_exposure_ms: Optional[float] = Field(description="Default exposure time")
    default_gain: Optional[float] = Field(description="Default gain")


# ============================================================================
# Image QC schemas
# ============================================================================


class ImageQCResponse(BaseSchema):
    """Response schema for single image QC."""
    
    image_id: UUID = Field(description="Image ID")
    focus_score: float = Field(description="Focus quality score (0-1)")
    focus_algorithm: str = Field(description="Focus algorithm used")
    is_focused: bool = Field(description="Whether image passes focus threshold")
    saturation_percent: float = Field(description="Percentage of saturated pixels")
    is_saturated: bool = Field(description="Whether image has excessive saturation")
    uniformity_score: float = Field(description="Illumination uniformity score (0-1)")
    vignetting_detected: bool = Field(description="Whether vignetting was detected")
    artifact_count: Optional[int] = Field(None, description="Number of artifacts detected")
    artifact_percent: Optional[float] = Field(None, description="Percentage area with artifacts")
    overall_score: float = Field(description="Overall QC score (0-100)")
    passed_qc: bool = Field(description="Whether image passes all QC checks")
    issues: List[str] = Field(description="List of QC issues found")
    timestamp: datetime = Field(description="QC analysis timestamp")


class PlateQCResponse(BaseSchema):
    """Response schema for plate-wide QC report."""
    
    plate_id: UUID = Field(description="Plate ID")
    total_images: int = Field(description="Total images analyzed")
    passed_count: int = Field(description="Images that passed QC")
    failed_count: int = Field(description="Images that failed QC")
    average_focus_score: float = Field(description="Average focus score across plate")
    focus_heatmap: Dict[str, float] = Field(description="Focus scores by well position")
    saturation_alerts: List[str] = Field(description="Wells with saturation issues")
    uniformity_issues: List[str] = Field(description="Wells with uniformity problems")
    recommendations: List[str] = Field(description="QC improvement recommendations")


# ============================================================================
# 5D Browser schemas
# ============================================================================


class BrowseQuery(BaseSchema):
    """Request schema for 5D data browser."""
    
    plate_id: Optional[UUID] = Field(None, description="Filter by plate ID")
    well_position: Optional[str] = Field(None, description="Filter by well position (e.g., A01)")
    channel: Optional[str] = Field(None, description="Filter by channel name")
    z_slice: Optional[int] = Field(None, ge=0, description="Filter by Z slice")
    timepoint: Optional[int] = Field(None, ge=0, description="Filter by timepoint")
    min_focus_score: Optional[float] = Field(None, ge=0, le=1, description="Minimum focus score")
    passed_qc_only: bool = Field(False, description="Only return QC-passed images")
    skip: int = Field(0, ge=0, description="Number of results to skip")
    limit: int = Field(100, ge=1, le=1000, description="Maximum results to return")


class ImageSummary(BaseSchema):
    """Summary schema for browse results."""
    
    id: UUID = Field(description="Image ID")
    well_position: Optional[str] = Field(description="Well position")
    channel: str = Field(description="Channel name")
    z_slice: int = Field(description="Z slice index")
    timepoint: int = Field(description="Timepoint index")
    width: int = Field(description="Image width in pixels")
    height: int = Field(description="Image height in pixels")
    thumbnail_url: Optional[str] = Field(None, description="Thumbnail image URL")
    focus_score: Optional[float] = Field(None, description="Focus quality score")
    passed_qc: Optional[bool] = Field(None, description="QC pass/fail status")
    acquired_at: Optional[datetime] = Field(None, description="Image acquisition time")


class BrowseResponse(BaseSchema):
    """Response schema for 5D browser query."""
    
    images: List[ImageSummary] = Field(description="Image results")
    total_count: int = Field(description="Total matching images")
    available_channels: List[str] = Field(description="Available channels in dataset")
    z_range: Tuple[int, int] = Field(description="Z slice range (min, max)")
    t_range: Tuple[int, int] = Field(description="Timepoint range (min, max)")
    plate_info: Optional[Dict[str, Any]] = Field(None, description="Plate metadata")


# ============================================================================
# Thumbnail Generation schemas
# ============================================================================


class ThumbnailRequest(BaseSchema):
    """Request schema for thumbnail generation."""
    
    image_id: UUID = Field(description="Image ID to generate thumbnail for")
    size: int = Field(256, ge=64, le=1024, description="Thumbnail size in pixels")
    contrast_enhancement: bool = Field(True, description="Apply contrast enhancement")
    format: str = Field("jpeg", description="Output format (jpeg, png)")


class ThumbnailResponse(BaseSchema):
    """Response schema for thumbnail generation."""
    
    image_id: UUID = Field(description="Source image ID")
    thumbnail_url: str = Field(description="URL to generated thumbnail")
    size: int = Field(description="Thumbnail size in pixels")
    format: str = Field(description="Image format")
    created_at: datetime = Field(description="Generation timestamp")


# ============================================================================
# Metadata Export schemas
# ============================================================================


class MetadataExportRequest(BaseSchema):
    """Request schema for metadata export."""
    
    plate_id: Optional[UUID] = Field(None, description="Export metadata for specific plate")
    fileset_id: Optional[UUID] = Field(None, description="Export metadata for import batch")
    format: str = Field("csv", description="Export format (csv, json, xlsx)")
    include_qc: bool = Field(True, description="Include QC metrics")
    include_acquisition: bool = Field(True, description="Include acquisition metadata")


class MetadataExportResponse(BaseSchema):
    """Response schema for metadata export."""
    
    export_id: UUID = Field(description="Export job ID")
    download_url: str = Field(description="URL to download exported file")
    format: str = Field(description="Export format")
    record_count: int = Field(description="Number of records exported")
    file_size_bytes: int = Field(description="Export file size")
    expires_at: datetime = Field(description="Download link expiration")


# Enable forward references
MicroscopeResponse.model_rebuild()
