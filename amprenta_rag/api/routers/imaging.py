"""Imaging analysis API endpoints."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status, UploadFile, File, Form
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select, func

from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.async_dependencies import get_async_database_session
from amprenta_rag.api import schemas as api_schemas
from amprenta_rag.api.schemas_imaging.imaging import (
    OMETiffImportRequest,
    OMETiffImportResponse,
    BatchImportRequest,
    BatchImportResponse,
    ImportStatusResponse,
    MicroscopeCreate,
    MicroscopeResponse,
    ObjectiveResponse,
    ChannelConfigResponse,
    ImageQCResponse,
    PlateQCResponse,
    BrowseQuery,
    BrowseResponse,
)
from amprenta_rag.models.auth import User
from amprenta_rag.models.chemistry import HTSWell
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.models_metadata import (
    Microscope, Objective, ChannelConfig, ImageFileSet
)
from amprenta_rag.imaging.cellpose_service import CellPoseService
from amprenta_rag.imaging.feature_extraction import FeatureExtractor
from amprenta_rag.imaging.storage import ImageStorage
from amprenta_rag.imaging.ome_parser import parse_ome_tiff
from amprenta_rag.imaging.vendor_parsers import parse_vendor_export
from amprenta_rag.imaging.image_qc import run_image_qc, generate_plate_qc_report
from amprenta_rag.jobs.tasks.imaging import process_batch_segmentation
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

router = APIRouter()

# Global service instances (would be configured from settings in production)
_cellpose_service: Optional[CellPoseService] = None
_feature_extractor: Optional[FeatureExtractor] = None
_image_storage: Optional[ImageStorage] = None


def get_cellpose_service() -> CellPoseService:
    """Get or create CellPose service instance."""
    global _cellpose_service
    if _cellpose_service is None:
        _cellpose_service = CellPoseService(
            model_type="cyto",
            gpu=True,
            tile_size=1024,
            overlap=128
        )
    return _cellpose_service


def get_feature_extractor() -> FeatureExtractor:
    """Get or create feature extractor instance."""
    global _feature_extractor
    if _feature_extractor is None:
        _feature_extractor = FeatureExtractor()
    return _feature_extractor


def get_image_storage() -> ImageStorage:
    """Get or create image storage instance."""
    global _image_storage
    if _image_storage is None:
        # In production, this would be configured from settings
        _image_storage = ImageStorage.create_local("data/imaging")
    return _image_storage


@router.post(
    "/imaging/upload",
    response_model=api_schemas.ImageUploadResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Upload microscopy image",
    description="Upload a microscopy image with metadata and optional well association."
)
async def upload_image(
    file: UploadFile = File(..., description="Image file (TIFF, PNG, JPG)"),
    well_id: Optional[str] = Form(None, description="Well ID to associate with image"),
    channel: str = Form(..., description="Channel name (e.g., DAPI, GFP, RFP)"),
    z_slice: int = Form(0, description="Z-stack slice index"),
    timepoint: int = Form(0, description="Time series index"),
    pixel_size_um: Optional[float] = Form(None, description="Pixel size in microns"),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    storage: ImageStorage = Depends(get_image_storage)
) -> ImageUploadResponse:
    """Upload microscopy image with metadata."""
    try:
        # Validate file type
        if file.content_type not in ["image/tiff", "image/png", "image/jpeg", "image/jpg"]:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Unsupported file type. Use TIFF, PNG, or JPG."
            )
        
        # Validate well_id if provided
        well_uuid = None
        if well_id:
            try:
                well_uuid = UUID(well_id)
                # Verify well exists
                from amprenta_rag.models.chemistry import HTSWell
                result = await db.execute(
                    select(HTSWell).filter(HTSWell.id == well_uuid)
                )
                well = result.scalars().first()
                if not well:
                    raise HTTPException(
                        status_code=status.HTTP_404_NOT_FOUND,
                        detail=f"Well {well_id} not found"
                    )
            except ValueError:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail="Invalid well_id format"
                )
        
        # Read image data
        image_data = file.file.read()
        
        # Save image to storage
        image_path = storage.save_image(
            image_data=image_data,
            well_id=well_id or "upload",
            channel=channel,
            z_slice=z_slice,
            timepoint=timepoint,
            format=file.content_type.split("/")[1]
        )
        
        # Get image metadata
        from PIL import Image
        import io
        
        pil_image = Image.open(io.BytesIO(image_data))
        width, height = pil_image.size
        
        # Create database record
        image_record = MicroscopyImage(
            well_id=well_uuid,
            channel=channel,
            z_slice=z_slice,
            timepoint=timepoint,
            width=width,
            height=height,
            bit_depth=8 if pil_image.mode == "L" else 16,
            pixel_size_um=pixel_size_um,
            image_path=image_path,
            image_metadata={
                "original_filename": file.filename,
                "content_type": file.content_type,
                "file_size": len(image_data)
            }
        )
        
        db.add(image_record)
        await db.commit()
        await db.refresh(image_record)
        
        # Generate initial thumbnail
        try:
            from amprenta_rag.imaging.thumbnails import get_thumbnail_service
            thumbnail_service = get_thumbnail_service()
            thumbnail_path = thumbnail_service.get_or_generate(image_record)
            if thumbnail_path:
                image_record.thumbnail_path = thumbnail_path
                await db.commit()
        except Exception as e:
            logger.warning(f"Failed to generate initial thumbnail: {e}")
            # Non-blocking - continue without thumbnail
        
        logger.info(f"Uploaded image {image_record.id} for user {current_user.id}")
        
        return api_schemas.ImageUploadResponse(
            image_id=image_record.id,
            image_path=image_path,
            width=width,
            height=height,
            channel=channel,
            well_id=well_uuid,
            message="Image uploaded successfully"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to upload image: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to upload image"
        )


@router.post(
    "/imaging/segment",
    response_model=api_schemas.SegmentationResponse,
    status_code=status.HTTP_200_OK,
    summary="Segment cells in image",
    description="Run CellPose segmentation on a microscopy image."
)
async def segment_image(
    request: SegmentationRequest,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    cellpose_service: CellPoseService = Depends(get_cellpose_service),
    feature_extractor: FeatureExtractor = Depends(get_feature_extractor),
    storage: ImageStorage = Depends(get_image_storage)
) -> SegmentationResponse:
    """Segment cells in microscopy image using CellPose."""
    try:
        # Get image record
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.id == request.image_id)
        )
        image = result.scalars().first()
        if not image:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Image {request.image_id} not found"
            )
        
        # Load image from storage
        image_array = storage.load_image(image.image_path)
        
        # Validate image for segmentation
        if not cellpose_service.validate_image(image_array):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Invalid image for segmentation"
            )
        
        # Run segmentation
        masks, flows = cellpose_service.segment(
            image=image_array,
            diameter=request.diameter,
            channels=request.channels
        )
        
        # Count cells
        cell_count = cellpose_service.count_cells(masks)
        
        # Save segmentation mask
        mask_path = storage.save_mask(
            mask_data=masks,
            segmentation_id=str(request.image_id)
        )
        
        # Create segmentation record
        segmentation = CellSegmentation(
            image_id=request.image_id,
            model_name=request.model_name or "cyto",
            model_version="3.0",
            cell_count=cell_count,
            mask_path=mask_path,
            parameters={
                "diameter": request.diameter,
                "channels": request.channels,
                "model_name": request.model_name
            },
            confidence_score=0.9  # Placeholder - CellPose doesn't return confidence
        )
        
        db.add(segmentation)
        await db.commit()
        await db.refresh(segmentation)
        
        # Extract features if requested
        features_extracted = False
        if request.extract_features:
            try:
                # Extract morphological features
                morphology_features = feature_extractor.extract_morphology_features(masks)
                
                # Save features to database
                for i, feature in enumerate(morphology_features):
                    cell_feature = CellFeature(
                        segmentation_id=segmentation.id,
                        cell_id=i + 1,
                        area=feature.area,
                        perimeter=feature.perimeter,
                        circularity=feature.circularity,
                        eccentricity=feature.eccentricity,
                        solidity=feature.solidity,
                        centroid_x=feature.centroid_x,
                        centroid_y=feature.centroid_y,
                        intensity_features={},  # No intensity features without multi-channel
                        texture_features={},
                        custom_features=feature.to_dict()
                    )
                    db.add(cell_feature)
                
                await db.commit()
                features_extracted = True
                
            except Exception as e:
                logger.warning(f"Failed to extract features: {e}")
        
        logger.info(f"Segmented image {request.image_id}, found {cell_count} cells")
        
        return api_schemas.SegmentationResponse(
            segmentation_id=segmentation.id,
            image_id=request.image_id,
            cell_count=cell_count,
            model_name=segmentation.model_name,
            mask_path=mask_path,
            features_extracted=features_extracted,
            processing_time_seconds=0.0  # Placeholder
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to segment image: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to segment image"
        )


@router.post(
    "/imaging/segment-batch",
    response_model=api_schemas.BatchSegmentationResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Queue batch segmentation",
    description="Queue multiple images for batch segmentation processing using Celery."
)
async def segment_batch(
    request: BatchSegmentationRequest,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> BatchSegmentationResponse:
    """Queue batch segmentation task."""
    try:
        # Validate that all images exist
        image_ids = [str(img_id) for img_id in request.image_ids]
        result = await db.execute(
            select(MicroscopyImage).filter(
                MicroscopyImage.id.in_(request.image_ids)
            )
        )
        existing_images = result.scalars().all()
        
        if len(existing_images) != len(request.image_ids):
            missing_ids = set(request.image_ids) - {img.id for img in existing_images}
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Images not found: {list(missing_ids)}"
            )
        
        # Queue Celery task
        task = process_batch_segmentation.delay(
            image_ids=image_ids,
            model_name=request.model_name,
            diameter=request.diameter,
            channels=request.channels,
            extract_features=request.extract_features,
            user_id=str(current_user.id)
        )
        
        logger.info(f"Queued batch segmentation task {task.id} for {len(image_ids)} images")
        
        return api_schemas.BatchSegmentationResponse(
            task_id=task.id,
            image_count=len(request.image_ids),
            status="queued",
            message=f"Batch segmentation queued for {len(request.image_ids)} images"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to queue batch segmentation: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to queue batch segmentation"
        )


@router.get(
    "/imaging/images/{image_id}",
    response_model=api_schemas.ImageMetadataResponse,
    summary="Get image metadata",
    description="Retrieve metadata for a microscopy image."
)
async def get_image_metadata(
    image_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> ImageMetadataResponse:
    """Get microscopy image metadata."""
    try:
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.id == image_id)
        )
        image = result.scalars().first()
        if not image:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Image {image_id} not found"
            )
        
        return api_schemas.ImageMetadataResponse(
            image_id=image.id,
            well_id=image.well_id,
            channel=image.channel,
            z_slice=image.z_slice,
            timepoint=image.timepoint,
            width=image.width,
            height=image.height,
            bit_depth=image.bit_depth,
            pixel_size_um=image.pixel_size_um,
            image_path=image.image_path,
            metadata=image.image_metadata,
            acquired_at=image.acquired_at,
            created_at=image.created_at
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get image metadata: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get image metadata"
        )


@router.get(
    "/imaging/segmentations/{segmentation_id}",
    response_model=api_schemas.SegmentationResultResponse,
    summary="Get segmentation results",
    description="Retrieve segmentation results and metadata."
)
async def get_segmentation(
    segmentation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> SegmentationResultResponse:
    """Get segmentation results."""
    try:
        result = await db.execute(
            select(CellSegmentation).filter(CellSegmentation.id == segmentation_id)
        )
        segmentation = result.scalars().first()
        if not segmentation:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Segmentation {segmentation_id} not found"
            )
        
        return api_schemas.SegmentationResultResponse(
            segmentation_id=segmentation.id,
            image_id=segmentation.image_id,
            model_name=segmentation.model_name,
            model_version=segmentation.model_version,
            cell_count=segmentation.cell_count,
            mask_path=segmentation.mask_path,
            parameters=segmentation.parameters,
            confidence_score=segmentation.confidence_score,
            created_at=segmentation.created_at
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get segmentation: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get segmentation"
        )


@router.get(
    "/imaging/features/{segmentation_id}",
    response_model=api_schemas.CellFeaturesResponse,
    summary="Get cell features",
    description="Retrieve extracted cell features for a segmentation."
)
async def get_features(
    segmentation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> CellFeaturesResponse:
    """Get cell features for segmentation."""
    try:
        # Check segmentation exists
        result = await db.execute(
            select(CellSegmentation).filter(CellSegmentation.id == segmentation_id)
        )
        segmentation = result.scalars().first()
        if not segmentation:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Segmentation {segmentation_id} not found"
            )
        
        # Get features
        result = await db.execute(
            select(CellFeature).filter(CellFeature.segmentation_id == segmentation_id)
        )
        features = result.scalars().all()
        
        # Convert to response format
        feature_list = []
        for feature in features:
            feature_dict = {
                "cell_id": feature.cell_id,
                "area": feature.area,
                "perimeter": feature.perimeter,
                "circularity": feature.circularity,
                "eccentricity": feature.eccentricity,
                "solidity": feature.solidity,
                "centroid_x": feature.centroid_x,
                "centroid_y": feature.centroid_y,
                "intensity_features": feature.intensity_features,
                "texture_features": feature.texture_features,
                "custom_features": feature.custom_features
            }
            feature_list.append(feature_dict)
        
        return api_schemas.CellFeaturesResponse(
            segmentation_id=segmentation_id,
            cell_count=len(features),
            features=feature_list
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get features: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get features"
        )


@router.get(
    "/imaging/wells/{well_id}/summary",
    response_model=api_schemas.WellSummaryResponse,
    summary="Get well summary",
    description="Get aggregated features and statistics for all images in a well."
)
async def get_well_summary(
    well_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    feature_extractor: FeatureExtractor = Depends(get_feature_extractor)
) -> WellSummaryResponse:
    """Get well-level aggregated imaging summary."""
    try:
        # Check well exists
        from amprenta_rag.models.chemistry import HTSWell
        result = await db.execute(
            select(HTSWell).filter(HTSWell.id == well_id)
        )
        well = result.scalars().first()
        if not well:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Well {well_id} not found"
            )
        
        # Get all images for this well
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.well_id == well_id)
        )
        images = result.scalars().all()
        
        if not images:
            return api_schemas.WellSummaryResponse(
                well_id=well_id,
                image_count=0,
                total_cell_count=0,
                channels=[],
                aggregated_features={},
                summary_metrics={}
            )
        
        # Get all segmentations for these images
        image_ids = [img.id for img in images]
        result = await db.execute(
            select(CellSegmentation).filter(CellSegmentation.image_id.in_(image_ids))
        )
        segmentations = result.scalars().all()
        
        # Get all features for these segmentations
        segmentation_ids = [seg.id for seg in segmentations]
        result = await db.execute(
            select(CellFeature).filter(CellFeature.segmentation_id.in_(segmentation_ids))
        )
        features = result.scalars().all()
        
        # Calculate summary statistics
        total_cell_count = sum(seg.cell_count for seg in segmentations)
        channels = list(set(img.channel for img in images))
        
        # Aggregate features if available
        aggregated_features = {}
        summary_metrics = {}
        
        if features:
            # Convert to feature extractor format for aggregation
            from amprenta_rag.imaging.feature_extraction import CellMorphologyFeatures
            
            morphology_features = []
            for feature in features:
                morph_feature = CellMorphologyFeatures(
                    area=feature.area or 0,
                    perimeter=feature.perimeter or 0,
                    major_axis_length=0,  # Not stored in our model
                    minor_axis_length=0,  # Not stored in our model
                    eccentricity=feature.eccentricity or 0,
                    solidity=feature.solidity or 0,
                    extent=0,  # Not stored
                    centroid_x=feature.centroid_x or 0,
                    centroid_y=feature.centroid_y or 0,
                    circularity=feature.circularity or 0,
                    aspect_ratio=1.0  # Default
                )
                morphology_features.append(morph_feature)
            
            # Aggregate to well level
            well_aggregated = feature_extractor.aggregate_to_well(morphology_features, [])
            aggregated_features = well_aggregated.to_dict()
            
            # Calculate high-level metrics
            summary_metrics = feature_extractor.calculate_well_metrics(well_aggregated)
        
        return api_schemas.WellSummaryResponse(
            well_id=well_id,
            image_count=len(images),
            total_cell_count=total_cell_count,
            channels=channels,
            aggregated_features=aggregated_features,
            summary_metrics=summary_metrics
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get well summary: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get well summary"
        )


# ============================================================================
# OME-TIFF Import Endpoints
# ============================================================================


@router.post(
    "/imaging/import/ome-tiff",
    response_model=OMETiffImportResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Import OME-TIFF file",
    description="Import OME-TIFF file with metadata extraction and optional well association."
)
async def import_ome_tiff(
    file: UploadFile = File(..., description="OME-TIFF file"),
    request: OMETiffImportRequest = Depends(),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    storage: ImageStorage = Depends(get_image_storage)
) -> OMETiffImportResponse:
    """Import OME-TIFF file with metadata extraction."""
    try:
        # Validate file type
        if not file.filename or not file.filename.lower().endswith(('.tif', '.tiff', '.ome.tif', '.ome.tiff')):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="File must be a TIFF or OME-TIFF file"
            )
        
        # Read file data
        file_data = file.file.read()
        
        # Save to temporary location for parsing
        import tempfile
        import os
        
        with tempfile.NamedTemporaryFile(suffix='.tiff', delete=False) as tmp_file:
            tmp_file.write(file_data)
            tmp_path = tmp_file.name
        
        try:
            # Parse OME-TIFF metadata
            ome_metadata = parse_ome_tiff(tmp_path)
            
            # Extract basic image info
            dimensions = {
                "x": ome_metadata.dimensions.size_x,
                "y": ome_metadata.dimensions.size_y,
                "z": ome_metadata.dimensions.size_z,
                "c": ome_metadata.dimensions.size_c,
                "t": ome_metadata.dimensions.size_t
            }
            
            # Get channel names
            channels = [ch.name or f"Channel_{i}" for i, ch in enumerate(ome_metadata.channels)]
            
            # Apply channel mapping if provided
            if request.channel_mapping:
                mapped_channels = []
                for i, original_name in enumerate(channels):
                    mapped_name = request.channel_mapping.get(str(i), original_name)
                    mapped_channels.append(mapped_name)
                channels = mapped_channels
            
            # Save image to storage
            image_path = storage.save_image(
                image_data=file_data,
                well_id=str(request.well_id) if request.well_id else "ome_import",
                channel=channels[0] if channels else "default",
                z_slice=0,
                timepoint=0,
                format="tiff"
            )
            
            # Create primary image record
            image_record = MicroscopyImage(
                well_id=request.well_id,
                channel=channels[0] if channels else "default",
                z_slice=0,
                timepoint=0,
                width=dimensions["x"],
                height=dimensions["y"],
                bit_depth=16,  # Default for OME-TIFF
                pixel_size_um=ome_metadata.dimensions.pixel_size_x_um,
                image_path=image_path,
                ome_metadata=ome_metadata.raw_json,
                image_metadata={
                    "original_filename": file.filename,
                    "ome_uuid": ome_metadata.ome_uuid,
                    "acquisition_date": ome_metadata.acquisition_date.isoformat() if ome_metadata.acquisition_date else None,
                    "description": ome_metadata.description
                }
            )
            
            # Set instrument info if available
            if ome_metadata.instrument:
                image_record.image_metadata["instrument"] = {
                    "microscope_name": ome_metadata.instrument.microscope_name,
                    "microscope_model": ome_metadata.instrument.microscope_model,
                    "objective_name": ome_metadata.instrument.objective_name,
                    "objective_magnification": ome_metadata.instrument.objective_magnification,
                    "objective_na": ome_metadata.instrument.objective_na
                }
            
            db.add(image_record)
            await db.commit()
            await db.refresh(image_record)
            
            logger.info(f"Imported OME-TIFF {file.filename} as image {image_record.id}")
            
            return OMETiffImportResponse(
                image_id=image_record.id,
                filename=file.filename,
                dimensions=dimensions,
                channels=channels,
                instrument=ome_metadata.instrument.microscope_name if ome_metadata.instrument else None,
                pixel_size_um=ome_metadata.dimensions.pixel_size_x_um,
                ome_metadata=ome_metadata.raw_json
            )
            
        finally:
            # Clean up temporary file
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to import OME-TIFF: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to import OME-TIFF: {str(e)}"
        )


@router.post(
    "/imaging/import/batch",
    response_model=BatchImportResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Batch import from vendor export",
    description="Import multiple images from vendor export directory (Opera, ImageXpress, Cell Voyager)."
)
async def import_batch_vendor(
    request: BatchImportRequest,
    import_path: str = Form(..., description="Path to vendor export directory"),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> BatchImportResponse:
    """Import batch of images from vendor export."""
    try:
        # Parse vendor export
        import_result = parse_vendor_export(import_path, vendor=request.vendor)
        
        if import_result.errors:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Failed to parse vendor export: {'; '.join(import_result.errors)}"
            )
        
        # Create ImageFileSet record
        fileset = ImageFileSet(
            plate_id=request.plate_id,
            vendor=import_result.vendor,
            import_path=import_path,
            file_count=len(import_result.images),
            image_count=len(import_result.images),
            import_status="pending"
        )
        
        db.add(fileset)
        await db.commit()
        await db.refresh(fileset)
        
        # NOTE: Background import jobs tracked in ROADMAP
        # For now, just update status to importing
        fileset.import_status = "importing"
        await db.commit()
        
        logger.info(f"Queued batch import {fileset.id} for {len(import_result.images)} images")
        
        return BatchImportResponse(
            fileset_id=fileset.id,
            status=fileset.import_status,
            vendor=import_result.vendor,
            total_images=len(import_result.images),
            imported_count=0,
            errors=import_result.errors
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to start batch import: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to start batch import: {str(e)}"
        )


@router.get(
    "/imaging/import/{fileset_id}/status",
    response_model=ImportStatusResponse,
    summary="Get import job status",
    description="Get status and progress of batch import job."
)
async def get_import_status(
    fileset_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> ImportStatusResponse:
    """Get import job status."""
    try:
        result = await db.execute(
            select(ImageFileSet).filter(ImageFileSet.id == fileset_id)
        )
        fileset = result.scalars().first()
        if not fileset:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Import job {fileset_id} not found"
            )
        
        # Calculate progress
        progress_percent = 0.0
        if fileset.file_count > 0:
            progress_percent = (fileset.image_count / fileset.file_count) * 100.0
        
        return ImportStatusResponse(
            fileset_id=fileset.id,
            status=fileset.import_status,
            progress_percent=progress_percent,
            total_images=fileset.file_count,
            imported_count=fileset.image_count,
            failed_count=0,  # NOTE: Failed import tracking tracked in ROADMAP
            errors=[fileset.error_message] if fileset.error_message else [],
            warnings=[],
            estimated_completion=None
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get import status: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get import status"
        )


# ============================================================================
# Instrument Management Endpoints
# ============================================================================


@router.get(
    "/imaging/instruments",
    response_model=List[MicroscopeResponse],
    summary="List microscopes",
    description="Get list of registered microscope instruments."
)
async def list_microscopes(
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> List[MicroscopeResponse]:
    """List all registered microscopes."""
    try:
        result = await db.execute(
            select(Microscope).filter(Microscope.is_active.is_(True))
        )
        microscopes = result.scalars().all()
        
        result = []
        for microscope in microscopes:
            # Get associated objectives
            obj_result = await db.execute(
                select(Objective).filter(
                    Objective.microscope_id == microscope.id,
                    Objective.is_active.is_(True)
                )
            )
            objectives = obj_result.scalars().all()
            
            # Get associated channel configs
            chan_result = await db.execute(
                select(ChannelConfig).filter(ChannelConfig.microscope_id == microscope.id)
            )
            channels = chan_result.scalars().all()
            
            result.append(MicroscopeResponse(
                id=microscope.id,
                name=microscope.name,
                manufacturer=microscope.manufacturer,
                model=microscope.model,
                serial_number=microscope.serial_number,
                facility_location=microscope.facility_location,
                is_active=microscope.is_active,
                objectives=[
                    ObjectiveResponse(
                        id=obj.id,
                        microscope_id=obj.microscope_id,
                        name=obj.name,
                        magnification=obj.magnification,
                        numerical_aperture=obj.numerical_aperture,
                        immersion=obj.immersion,
                        working_distance_mm=obj.working_distance_mm,
                        correction=obj.correction,
                        is_active=obj.is_active
                    ) for obj in objectives
                ],
                channels=[
                    ChannelConfigResponse(
                        id=ch.id,
                        microscope_id=ch.microscope_id,
                        channel_name=ch.channel_name,
                        fluorophore=ch.fluorophore,
                        excitation_wavelength_nm=None,  # NOTE: Filter set wavelengths tracked in ROADMAP
                        emission_wavelength_nm=None,    # NOTE: Filter set wavelengths tracked in ROADMAP
                        default_exposure_ms=ch.default_exposure_ms,
                        default_gain=ch.default_gain
                    ) for ch in channels
                ],
                created_at=microscope.created_at
            ))
        
        return result
        
    except Exception as e:
        logger.error(f"Failed to list microscopes: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to list microscopes"
        )


@router.post(
    "/imaging/instruments",
    response_model=MicroscopeResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Register microscope",
    description="Register a new microscope instrument."
)
async def create_microscope(
    request: MicroscopeCreate,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> MicroscopeResponse:
    """Register new microscope."""
    try:
        # Check for duplicate serial number
        if request.serial_number:
            result = await db.execute(
                select(Microscope).filter(Microscope.serial_number == request.serial_number)
            )
            existing = result.scalars().first()
            if existing:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"Microscope with serial number {request.serial_number} already exists"
                )
        
        microscope = Microscope(
            name=request.name,
            manufacturer=request.manufacturer,
            model=request.model,
            serial_number=request.serial_number,
            facility_location=request.facility_location
        )
        
        db.add(microscope)
        await db.commit()
        await db.refresh(microscope)
        
        logger.info(f"Created microscope {microscope.id}: {microscope.name}")
        
        return MicroscopeResponse(
            id=microscope.id,
            name=microscope.name,
            manufacturer=microscope.manufacturer,
            model=microscope.model,
            serial_number=microscope.serial_number,
            facility_location=microscope.facility_location,
            is_active=microscope.is_active,
            objectives=[],
            channels=[],
            created_at=microscope.created_at
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to create microscope: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to create microscope"
        )


@router.get(
    "/imaging/instruments/{microscope_id}",
    response_model=MicroscopeResponse,
    summary="Get microscope details",
    description="Get detailed information about a specific microscope."
)
async def get_microscope(
    microscope_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> MicroscopeResponse:
    """Get microscope details."""
    try:
        result = await db.execute(
            select(Microscope).filter(Microscope.id == microscope_id)
        )
        microscope = result.scalars().first()
        if not microscope:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Microscope {microscope_id} not found"
            )
        
        # Get objectives and channels (same as list_microscopes)
        obj_result = await db.execute(
            select(Objective).filter(
                Objective.microscope_id == microscope.id,
                Objective.is_active.is_(True)
            )
        )
        objectives = obj_result.scalars().all()
        
        chan_result = await db.execute(
            select(ChannelConfig).filter(ChannelConfig.microscope_id == microscope.id)
        )
        channels = chan_result.scalars().all()
        
        return MicroscopeResponse(
            id=microscope.id,
            name=microscope.name,
            manufacturer=microscope.manufacturer,
            model=microscope.model,
            serial_number=microscope.serial_number,
            facility_location=microscope.facility_location,
            is_active=microscope.is_active,
            objectives=[
                ObjectiveResponse(
                    id=obj.id,
                    microscope_id=obj.microscope_id,
                    name=obj.name,
                    magnification=obj.magnification,
                    numerical_aperture=obj.numerical_aperture,
                    immersion=obj.immersion,
                    working_distance_mm=obj.working_distance_mm,
                    correction=obj.correction,
                    is_active=obj.is_active
                ) for obj in objectives
            ],
            channels=[
                ChannelConfigResponse(
                    id=ch.id,
                    microscope_id=ch.microscope_id,
                    channel_name=ch.channel_name,
                    fluorophore=ch.fluorophore,
                    excitation_wavelength_nm=None,
                    emission_wavelength_nm=None,
                    default_exposure_ms=ch.default_exposure_ms,
                    default_gain=ch.default_gain
                ) for ch in channels
            ],
            created_at=microscope.created_at
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get microscope: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get microscope"
        )


@router.get(
    "/imaging/objectives",
    response_model=List[ObjectiveResponse],
    summary="List objectives",
    description="Get list of available objective lenses."
)
async def list_objectives(
    microscope_id: Optional[UUID] = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> List[ObjectiveResponse]:
    """List objectives, optionally filtered by microscope."""
    try:
        query = select(Objective).filter(Objective.is_active.is_(True))
        
        if microscope_id:
            query = query.filter(Objective.microscope_id == microscope_id)
        
        result = await db.execute(query)
        objectives = result.scalars().all()
        
        return [
            ObjectiveResponse(
                id=obj.id,
                microscope_id=obj.microscope_id,
                name=obj.name,
                magnification=obj.magnification,
                numerical_aperture=obj.numerical_aperture,
                immersion=obj.immersion,
                working_distance_mm=obj.working_distance_mm,
                correction=obj.correction,
                is_active=obj.is_active
            ) for obj in objectives
        ]
        
    except Exception as e:
        logger.error(f"Failed to list objectives: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to list objectives"
        )


@router.get(
    "/imaging/channels",
    response_model=List[ChannelConfigResponse],
    summary="List channel configurations",
    description="Get list of channel configurations for microscopes."
)
async def list_channel_configs(
    microscope_id: Optional[UUID] = None,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> List[ChannelConfigResponse]:
    """List channel configurations."""
    try:
        query = select(ChannelConfig)
        
        if microscope_id:
            query = query.filter(ChannelConfig.microscope_id == microscope_id)
        
        result = await db.execute(query)
        channels = result.scalars().all()
        
        return [
            ChannelConfigResponse(
                id=ch.id,
                microscope_id=ch.microscope_id,
                channel_name=ch.channel_name,
                fluorophore=ch.fluorophore,
                excitation_wavelength_nm=None,  # NOTE: Filter set wavelengths tracked in ROADMAP
                emission_wavelength_nm=None,    # NOTE: Filter set wavelengths tracked in ROADMAP
                default_exposure_ms=ch.default_exposure_ms,
                default_gain=ch.default_gain
            ) for ch in channels
        ]
        
    except Exception as e:
        logger.error(f"Failed to list channel configs: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to list channel configs"
        )


# ============================================================================
# Image QC Endpoints
# ============================================================================


@router.get(
    "/imaging/qc/image/{image_id}",
    response_model=ImageQCResponse,
    summary="Get single image QC",
    description="Get quality control metrics for a single image."
)
async def get_image_qc(
    image_id: UUID,
    run_artifact_detection: bool = False,
    persist: bool = False,  # NEW: Save results to DB
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    storage: ImageStorage = Depends(get_image_storage)
) -> ImageQCResponse:
    """Get QC metrics for single image."""
    try:
        # Get image record
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.id == image_id)
        )
        image = result.scalars().first()
        if not image:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Image {image_id} not found"
            )
        
        # Load image data
        image_array = storage.load_image(image.image_path)
        
        # Run QC analysis
        qc_result = run_image_qc(
            image_array,
            image_path=image.image_path,
            run_artifact_detection=run_artifact_detection
        )
        
        # Optionally persist QC results
        if persist:
            from amprenta_rag.imaging.models import ImageQCRecord
            qc_record = ImageQCRecord(
                image_id=image_id,
                focus_score=qc_result.focus.score,
                focus_algorithm=qc_result.focus.algorithm,
                is_focused=qc_result.focus.is_focused,
                saturation_percent=qc_result.saturation.saturated_percent,
                is_saturated=qc_result.saturation.is_saturated,
                uniformity_score=qc_result.uniformity.uniformity_score,
                vignetting_detected=qc_result.uniformity.vignetting_detected,
                artifact_count=qc_result.artifacts.artifact_count if qc_result.artifacts else None,
                artifact_percent=qc_result.artifacts.artifact_percent if qc_result.artifacts else None,
                overall_score=qc_result.overall_score,
                passed_qc=qc_result.passed_qc,
                issues=qc_result.issues,
            )
            db.add(qc_record)
            await db.commit()

        return ImageQCResponse(
            image_id=image_id,
            focus_score=qc_result.focus.score,
            focus_algorithm=qc_result.focus.algorithm,
            is_focused=qc_result.focus.is_focused,
            saturation_percent=qc_result.saturation.saturated_percent,
            is_saturated=qc_result.saturation.is_saturated,
            uniformity_score=qc_result.uniformity.uniformity_score,
            vignetting_detected=qc_result.uniformity.vignetting_detected,
            artifact_count=qc_result.artifacts.artifact_count if qc_result.artifacts else None,
            artifact_percent=qc_result.artifacts.artifact_percent if qc_result.artifacts else None,
            overall_score=qc_result.overall_score,
            passed_qc=qc_result.passed_qc,
            issues=qc_result.issues,
            timestamp=qc_result.timestamp
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get image QC: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get image QC"
        )


@router.get(
    "/imaging/qc/history/{image_id}",
    response_model=List[ImageQCResponse],
    summary="Get QC history for image",
    description="Get historical QC records for a microscopy image."
)
async def get_qc_history(
    image_id: UUID,
    limit: int = Query(10, ge=1, le=100, description="Maximum records to return"),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> List[ImageQCResponse]:
    """Get historical QC records for an image."""
    try:
        from amprenta_rag.imaging.models import ImageQCRecord
        
        result = await db.execute(
            select(ImageQCRecord)
            .filter(ImageQCRecord.image_id == image_id)
            .order_by(ImageQCRecord.created_at.desc())
            .limit(limit)
        )
        records = result.scalars().all()
        
        return [
            ImageQCResponse(
                image_id=record.image_id,
                focus_score=record.focus_score,
                focus_algorithm=record.focus_algorithm,
                is_focused=record.is_focused,
                saturation_percent=record.saturation_percent,
                is_saturated=record.is_saturated,
                uniformity_score=record.uniformity_score,
                vignetting_detected=record.vignetting_detected,
                artifact_count=record.artifact_count,
                artifact_percent=record.artifact_percent,
                overall_score=record.overall_score,
                passed_qc=record.passed_qc,
                issues=record.issues or [],
                timestamp=record.created_at
            )
            for record in records
        ]
        
    except Exception as e:
        logger.error(f"Failed to get QC history: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get QC history"
        )


@router.get(
    "/imaging/qc/plate/{plate_id}",
    response_model=PlateQCResponse,
    summary="Get plate QC report",
    description="Get quality control report for an entire plate."
)
async def get_plate_qc_report(
    plate_id: UUID,
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session),
    storage: ImageStorage = Depends(get_image_storage)
) -> PlateQCResponse:
    """Get plate-wide QC report."""
    try:
        # Get all images for the plate
        well_subquery = select(HTSWell.id).filter(HTSWell.plate_id == plate_id)
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.well_id.in_(well_subquery))
        )
        images = result.scalars().all()
        
        if not images:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"No images found for plate {plate_id}"
            )
        
        # Load images and generate report
        image_data = []
        for image in images:
            try:
                image_array = storage.load_image(image.image_path)
                # Resolve well position for images
                if image.well_id:
                    well_result = await db.execute(
                        select(HTSWell).filter(HTSWell.id == image.well_id)
                    )
                    well = well_result.scalars().first()
                    well_position = well.position if well else f"orphan_{image.id.hex[:8]}"
                else:
                    well_position = f"standalone_{image.id.hex[:8]}"
                image_data.append((well_position, image_array))
            except Exception as e:
                logger.warning(f"Failed to load image {image.id}: {e}")
                continue
        
        if not image_data:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to load any images for QC analysis"
            )
        
        # Generate plate QC report
        qc_report = generate_plate_qc_report(image_data, str(plate_id))
        
        return PlateQCResponse(
            plate_id=plate_id,
            total_images=qc_report.total_images,
            passed_count=qc_report.passed_count,
            failed_count=qc_report.failed_count,
            average_focus_score=qc_report.average_focus_score,
            focus_heatmap=qc_report.focus_heatmap,
            saturation_alerts=qc_report.saturation_alerts,
            uniformity_issues=qc_report.uniformity_issues,
            recommendations=qc_report.recommendations
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get plate QC report: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get plate QC report"
        )


# ============================================================================
# 5D Browser Endpoint
# ============================================================================


@router.get(
    "/imaging/thumbnails/{image_id}",
    summary="Get image thumbnail",
    description="Serve cached thumbnail for a microscopy image."
)
async def get_thumbnail(
    image_id: UUID,
    size: int = Query(256, ge=64, le=512, description="Thumbnail size in pixels"),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
):
    """Serve cached or generated thumbnail."""
    from fastapi.responses import FileResponse
    from amprenta_rag.imaging.thumbnails import get_thumbnail_service
    
    try:
        # Verify image exists
        result = await db.execute(
            select(MicroscopyImage).filter(MicroscopyImage.id == image_id)
        )
        image = result.scalars().first()
        if not image:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Image {image_id} not found"
            )
        
        # Get or generate thumbnail
        thumbnail_service = get_thumbnail_service()
        thumbnail_path = thumbnail_service.get_or_generate(image, size=size)
        
        if not thumbnail_path:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to generate thumbnail"
            )
        
        return FileResponse(
            thumbnail_path,
            media_type="image/jpeg",
            headers={"Cache-Control": "public, max-age=86400"}  # 24h cache
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get thumbnail: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to get thumbnail"
        )


@router.get(
    "/imaging/browse",
    response_model=BrowseResponse,
    summary="5D data browser query",
    description="Browse and filter microscopy images across 5 dimensions (X, Y, Z, C, T)."
)
async def browse_images(
    query: BrowseQuery = Depends(),
    current_user: User = Depends(get_current_user),
    db: AsyncSession = Depends(get_async_database_session)
) -> BrowseResponse:
    """Browse images with 5D filtering."""
    try:
        # Build query
        image_query = select(MicroscopyImage)
        
        # Apply filters
        if query.plate_id:
            well_subquery = select(HTSWell.id).filter(HTSWell.plate_id == query.plate_id)
            image_query = image_query.filter(MicroscopyImage.well_id.in_(well_subquery))
        
        if query.well_position:
            # NOTE: Well position filtering tracked in ROADMAP
            pass
        
        if query.channel:
            image_query = image_query.filter(MicroscopyImage.channel == query.channel)
        
        if query.z_slice is not None:
            image_query = image_query.filter(MicroscopyImage.z_slice == query.z_slice)
        
        if query.timepoint is not None:
            image_query = image_query.filter(MicroscopyImage.timepoint == query.timepoint)
        
        # Get total count before pagination
        count_result = await db.execute(select(func.count()).select_from(image_query.subquery()))
        total_count = count_result.scalar()
        
        # Apply pagination
        paginated_query = image_query.offset(query.skip).limit(query.limit)
        result = await db.execute(paginated_query)
        images = result.scalars().all()
        
        # Get available channels and ranges
        all_images_query = select(MicroscopyImage)
        if query.plate_id:
            well_subquery = select(HTSWell.id).filter(HTSWell.plate_id == query.plate_id)
            all_images_query = all_images_query.filter(MicroscopyImage.well_id.in_(well_subquery))
        
        # Get distinct channels
        channels_result = await db.execute(
            select(MicroscopyImage.channel.distinct())
        )
        available_channels = [row[0] for row in channels_result.all()]
        
        # Get z range
        z_range_result = await db.execute(
            select(func.min(MicroscopyImage.z_slice), func.max(MicroscopyImage.z_slice))
        )
        z_range_row = z_range_result.first()
        z_range = (z_range_row[0] or 0, z_range_row[1] or 0)
        
        # Get timepoint range
        t_range_result = await db.execute(
            select(func.min(MicroscopyImage.timepoint), func.max(MicroscopyImage.timepoint))
        )
        t_range_row = t_range_result.first()
        t_range = (t_range_row[0] or 0, t_range_row[1] or 0)
        
        # Convert to response format
        image_summaries = []
        for img in images:
            # Resolve well position for images
            if img.well_id:
                well_result = await db.execute(
                    select(HTSWell).filter(HTSWell.id == img.well_id)
                )
                well = well_result.scalars().first()
                well_position = well.position if well else f"orphan_{img.id.hex[:8]}"
            else:
                well_position = f"standalone_{img.id.hex[:8]}"
            
            image_summaries.append(ImageSummary(
                id=img.id,
                well_position=well_position,
                channel=img.channel,
                z_slice=img.z_slice,
                timepoint=img.timepoint,
                width=img.width,
                height=img.height,
                thumbnail_url=f"/imaging/thumbnails/{img.id}",
                focus_score=None,    # NOTE: QC data integration tracked in ROADMAP
                passed_qc=None,      # NOTE: QC data integration tracked in ROADMAP
                acquired_at=img.acquired_at
            ))
        
        return BrowseResponse(
            images=image_summaries,
            total_count=total_count,
            available_channels=available_channels,
            z_range=z_range,
            t_range=t_range,
            plate_info=None  # NOTE: Plate metadata tracked in ROADMAP
        )
        
    except Exception as e:
        logger.error(f"Failed to browse images: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to browse images"
        )
