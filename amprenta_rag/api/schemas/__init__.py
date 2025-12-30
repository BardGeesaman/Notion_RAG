"""API schemas package."""

# Simple approach: Define the essential schemas here to avoid circular imports
from datetime import datetime
from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field


class BaseSchema(BaseModel):
    """Base schema with common configuration."""

    model_config = ConfigDict(
        from_attributes=True,
        json_encoders={
            UUID: str,
            datetime: lambda v: v.isoformat(),
        }
    )


class AnnotationCreate(BaseSchema):
    """Schema for creating an annotation/note on an entity."""

    text: str
    annotation_type: Optional[str] = None


class ErrorResponse(BaseSchema):
    """Schema for error responses."""
    
    error: str = Field(..., description="Error message")
    detail: Optional[str] = Field(None, description="Detailed error information")


# Imaging schemas (existing ones needed for backward compatibility)
class ImageUploadRequest(BaseSchema):
    """Request for uploading microscopy image."""
    
    well_id: Optional[UUID] = Field(None, description="Well ID to associate with image")
    channel: str = Field(..., description="Channel name (e.g., DAPI, GFP, RFP)")
    z_slice: int = Field(0, description="Z-stack slice index")
    timepoint: int = Field(0, description="Time series index")
    pixel_size_um: Optional[float] = Field(None, description="Pixel size in microns")


class ImageUploadResponse(BaseSchema):
    """Response for image upload."""
    
    image_id: UUID = Field(description="Unique image identifier")
    image_path: str = Field(description="Storage path of uploaded image")
    width: int = Field(description="Image width in pixels")
    height: int = Field(description="Image height in pixels")
    channel: str = Field(description="Channel name")
    well_id: Optional[UUID] = Field(description="Associated well ID")
    message: str = Field(description="Success message")


class SegmentationRequest(BaseSchema):
    """Request for cell segmentation."""
    
    image_id: UUID = Field(description="Image to segment")
    model_name: Optional[str] = Field("cyto", description="CellPose model type")
    diameter: Optional[float] = Field(30.0, description="Expected cell diameter in pixels")
    channels: Optional[List[int]] = Field([0, 0], description="Channel configuration [cytoplasm, nucleus]")
    extract_features: bool = Field(True, description="Whether to extract morphological features")


class SegmentationResponse(BaseSchema):
    """Response for cell segmentation."""
    
    segmentation_id: UUID = Field(description="Unique segmentation identifier")
    image_id: UUID = Field(description="Source image ID")
    cell_count: int = Field(description="Number of cells detected")
    model_name: str = Field(description="Model used for segmentation")
    mask_path: str = Field(description="Storage path of segmentation mask")
    features_extracted: bool = Field(description="Whether features were extracted")
    processing_time_seconds: float = Field(description="Processing time")


class BatchSegmentationRequest(BaseSchema):
    """Request for batch segmentation."""
    
    image_ids: List[UUID] = Field(description="List of images to segment")
    model_name: Optional[str] = Field("cyto", description="CellPose model type")
    diameter: Optional[float] = Field(30.0, description="Expected cell diameter in pixels")
    channels: Optional[List[int]] = Field([0, 0], description="Channel configuration")
    extract_features: bool = Field(True, description="Whether to extract features")


class BatchSegmentationResponse(BaseSchema):
    """Response for batch segmentation."""
    
    task_id: str = Field(description="Celery task ID for tracking")
    image_count: int = Field(description="Number of images queued")
    status: str = Field(description="Task status")
    message: str = Field(description="Status message")


class ImageMetadataResponse(BaseSchema):
    """Response for image metadata."""
    
    image_id: UUID = Field(description="Image identifier")
    well_id: Optional[UUID] = Field(description="Associated well ID")
    channel: str = Field(description="Channel name")
    z_slice: int = Field(description="Z-stack slice index")
    timepoint: int = Field(description="Time series index")
    width: int = Field(description="Image width in pixels")
    height: int = Field(description="Image height in pixels")
    bit_depth: int = Field(description="Image bit depth")
    pixel_size_um: Optional[float] = Field(description="Pixel size in microns")
    image_path: str = Field(description="Storage path")
    metadata: Optional[Dict[str, Any]] = Field(description="Additional metadata")
    acquired_at: Optional[datetime] = Field(description="Image acquisition time")
    created_at: datetime = Field(description="Record creation time")


class SegmentationResultResponse(BaseSchema):
    """Response for segmentation results."""
    
    segmentation_id: UUID = Field(description="Segmentation identifier")
    image_id: UUID = Field(description="Source image ID")
    model_name: str = Field(description="Model used")
    model_version: Optional[str] = Field(description="Model version")
    cell_count: int = Field(description="Number of cells detected")
    mask_path: str = Field(description="Segmentation mask path")
    parameters: Optional[Dict[str, Any]] = Field(description="Segmentation parameters")
    confidence_score: Optional[float] = Field(description="Overall confidence score")
    created_at: datetime = Field(description="Segmentation time")


class CellFeaturesResponse(BaseSchema):
    """Response for cell features."""
    
    segmentation_id: UUID = Field(description="Source segmentation ID")
    cell_count: int = Field(description="Number of cells")
    features: List[Dict[str, Any]] = Field(description="List of cell features")


class WellSummaryResponse(BaseSchema):
    """Response for well-level summary."""
    
    well_id: UUID = Field(description="Well identifier")
    image_count: int = Field(description="Number of images in well")
    total_cell_count: int = Field(description="Total cells across all images")
    channels: List[str] = Field(description="Available channels")
    aggregated_features: Dict[str, Any] = Field(description="Aggregated morphology features")
    summary_metrics: Dict[str, Any] = Field(description="High-level well metrics")


# Dataset schemas (needed by datasets router)
class DatasetCreate(BaseSchema):
    """Schema for creating a dataset."""
    name: str
    description: Optional[str] = None


class DatasetUpdate(BaseSchema):
    """Schema for updating a dataset."""
    name: Optional[str] = None
    description: Optional[str] = None


class Dataset(BaseSchema):
    """Dataset response schema."""
    id: UUID
    name: str
    description: Optional[str] = None
    created_at: datetime
    updated_at: datetime


# Dataset finder schemas
class DatasetFinderRequest(BaseSchema):
    """Request schema for AI dataset finder."""
    query: str = Field(..., min_length=1, max_length=500)


class DatasetFinderResponse(BaseSchema):
    """Response schema for AI dataset finder."""
    query: str
    results: List[Dict[str, Any]] = Field(default_factory=list)


class EnrichmentStatusResponse(BaseSchema):
    """Response schema for enrichment status."""
    enriched: bool
    available_fields: List[str] = Field(default_factory=list)


class MetadataEnrichmentResponse(BaseSchema):
    """Response schema for metadata enrichment."""
    dataset_id: UUID
    success: bool
    enriched_fields: List[str] = Field(default_factory=list)
