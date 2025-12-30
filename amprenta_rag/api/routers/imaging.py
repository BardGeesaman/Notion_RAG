"""Imaging analysis API endpoints."""

from __future__ import annotations

import logging
from typing import List, Optional, Dict, Any
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File, Form
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_db
from amprenta_rag.api.schemas import (
    ImageUploadRequest,
    ImageUploadResponse,
    SegmentationRequest,
    SegmentationResponse,
    BatchSegmentationRequest,
    BatchSegmentationResponse,
    ImageMetadataResponse,
    SegmentationResultResponse,
    CellFeaturesResponse,
    WellSummaryResponse,
    ErrorResponse,
)
from amprenta_rag.models.auth import User
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.cellpose_service import CellPoseService
from amprenta_rag.imaging.feature_extraction import FeatureExtractor
from amprenta_rag.imaging.storage import ImageStorage
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
    response_model=ImageUploadResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Upload microscopy image",
    description="Upload a microscopy image with metadata and optional well association."
)
def upload_image(
    file: UploadFile = File(..., description="Image file (TIFF, PNG, JPG)"),
    well_id: Optional[str] = Form(None, description="Well ID to associate with image"),
    channel: str = Form(..., description="Channel name (e.g., DAPI, GFP, RFP)"),
    z_slice: int = Form(0, description="Z-stack slice index"),
    timepoint: int = Form(0, description="Time series index"),
    pixel_size_um: Optional[float] = Form(None, description="Pixel size in microns"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
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
                well = db.query(HTSWell).filter(HTSWell.id == well_uuid).first()
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
        import numpy as np
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
        db.commit()
        db.refresh(image_record)
        
        logger.info(f"Uploaded image {image_record.id} for user {current_user.id}")
        
        return ImageUploadResponse(
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
    response_model=SegmentationResponse,
    status_code=status.HTTP_200_OK,
    summary="Segment cells in image",
    description="Run CellPose segmentation on a microscopy image."
)
def segment_image(
    request: SegmentationRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
    cellpose_service: CellPoseService = Depends(get_cellpose_service),
    feature_extractor: FeatureExtractor = Depends(get_feature_extractor),
    storage: ImageStorage = Depends(get_image_storage)
) -> SegmentationResponse:
    """Segment cells in microscopy image using CellPose."""
    try:
        # Get image record
        image = db.query(MicroscopyImage).filter(MicroscopyImage.id == request.image_id).first()
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
        db.commit()
        db.refresh(segmentation)
        
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
                
                db.commit()
                features_extracted = True
                
            except Exception as e:
                logger.warning(f"Failed to extract features: {e}")
        
        logger.info(f"Segmented image {request.image_id}, found {cell_count} cells")
        
        return SegmentationResponse(
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
    response_model=BatchSegmentationResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Queue batch segmentation",
    description="Queue multiple images for batch segmentation processing using Celery."
)
def segment_batch(
    request: BatchSegmentationRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> BatchSegmentationResponse:
    """Queue batch segmentation task."""
    try:
        # Validate that all images exist
        image_ids = [str(img_id) for img_id in request.image_ids]
        existing_images = db.query(MicroscopyImage).filter(
            MicroscopyImage.id.in_(request.image_ids)
        ).all()
        
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
        
        return BatchSegmentationResponse(
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
    response_model=ImageMetadataResponse,
    summary="Get image metadata",
    description="Retrieve metadata for a microscopy image."
)
def get_image_metadata(
    image_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> ImageMetadataResponse:
    """Get microscopy image metadata."""
    try:
        image = db.query(MicroscopyImage).filter(MicroscopyImage.id == image_id).first()
        if not image:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Image {image_id} not found"
            )
        
        return ImageMetadataResponse(
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
    response_model=SegmentationResultResponse,
    summary="Get segmentation results",
    description="Retrieve segmentation results and metadata."
)
def get_segmentation(
    segmentation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> SegmentationResultResponse:
    """Get segmentation results."""
    try:
        segmentation = db.query(CellSegmentation).filter(
            CellSegmentation.id == segmentation_id
        ).first()
        if not segmentation:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Segmentation {segmentation_id} not found"
            )
        
        return SegmentationResultResponse(
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
    response_model=CellFeaturesResponse,
    summary="Get cell features",
    description="Retrieve extracted cell features for a segmentation."
)
def get_features(
    segmentation_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> CellFeaturesResponse:
    """Get cell features for segmentation."""
    try:
        # Check segmentation exists
        segmentation = db.query(CellSegmentation).filter(
            CellSegmentation.id == segmentation_id
        ).first()
        if not segmentation:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Segmentation {segmentation_id} not found"
            )
        
        # Get features
        features = db.query(CellFeature).filter(
            CellFeature.segmentation_id == segmentation_id
        ).all()
        
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
        
        return CellFeaturesResponse(
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
    response_model=WellSummaryResponse,
    summary="Get well summary",
    description="Get aggregated features and statistics for all images in a well."
)
def get_well_summary(
    well_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
    feature_extractor: FeatureExtractor = Depends(get_feature_extractor)
) -> WellSummaryResponse:
    """Get well-level aggregated imaging summary."""
    try:
        # Check well exists
        from amprenta_rag.models.chemistry import HTSWell
        well = db.query(HTSWell).filter(HTSWell.id == well_id).first()
        if not well:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Well {well_id} not found"
            )
        
        # Get all images for this well
        images = db.query(MicroscopyImage).filter(
            MicroscopyImage.well_id == well_id
        ).all()
        
        if not images:
            return WellSummaryResponse(
                well_id=well_id,
                image_count=0,
                total_cell_count=0,
                channels=[],
                aggregated_features={},
                summary_metrics={}
            )
        
        # Get all segmentations for these images
        image_ids = [img.id for img in images]
        segmentations = db.query(CellSegmentation).filter(
            CellSegmentation.image_id.in_(image_ids)
        ).all()
        
        # Get all features for these segmentations
        segmentation_ids = [seg.id for seg in segmentations]
        features = db.query(CellFeature).filter(
            CellFeature.segmentation_id.in_(segmentation_ids)
        ).all()
        
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
        
        return WellSummaryResponse(
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
