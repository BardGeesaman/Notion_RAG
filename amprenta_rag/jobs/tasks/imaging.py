"""Celery tasks for imaging and microscopy processing."""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import List, Optional, Dict, Any
from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@celery_app.task(bind=True, max_retries=3, default_retry_delay=60, queue='default')
def process_batch_segmentation(
    self,
    image_ids: List[str],
    model_name: str = "cyto",
    diameter: Optional[float] = 30.0,
    channels: Optional[List[int]] = None,
    extract_features: bool = True,
    user_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Process batch cell segmentation for multiple images.
    
    Args:
        image_ids: List of image UUIDs to process
        model_name: CellPose model type
        diameter: Expected cell diameter in pixels
        channels: Channel configuration [cytoplasm, nucleus]
        extract_features: Whether to extract morphological features
        user_id: User who initiated the task
    
    Returns:
        Dictionary with processing results
    """
    if channels is None:
        channels = [0, 0]
    
    try:
        from amprenta_rag.database.session import db_session
        from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
        from amprenta_rag.imaging.cellpose_service import CellPoseService
        from amprenta_rag.imaging.feature_extraction import FeatureExtractor
        from amprenta_rag.imaging.storage import ImageStorage
        
        logger.info(f"Starting batch segmentation for {len(image_ids)} images")
        
        # Initialize services
        cellpose_service = CellPoseService(
            model_type=model_name,
            gpu=True,
            tile_size=1024,
            overlap=128
        )
        feature_extractor = FeatureExtractor()
        storage = ImageStorage.create_local("data/imaging")
        
        # Track results
        results = {
            "task_id": self.request.id,
            "total_images": len(image_ids),
            "processed_images": 0,
            "successful_segmentations": 0,
            "failed_segmentations": 0,
            "total_cells_found": 0,
            "features_extracted": 0,
            "errors": [],
            "processing_time_seconds": 0.0,
            "started_at": datetime.now(timezone.utc).isoformat(),
            "user_id": user_id
        }
        
        start_time = datetime.now()
        
        with db_session() as db:
            # Convert string IDs to UUIDs
            image_uuids = [UUID(img_id) for img_id in image_ids]
            
            # Get all images
            images = db.query(MicroscopyImage).filter(
                MicroscopyImage.id.in_(image_uuids)
            ).all()
            
            if len(images) != len(image_ids):
                found_ids = {str(img.id) for img in images}
                missing_ids = set(image_ids) - found_ids
                error_msg = f"Images not found: {list(missing_ids)}"
                logger.error(error_msg)
                results["errors"].append(error_msg)
                return results
            
            # Process each image
            for i, image in enumerate(images):
                try:
                    logger.info(f"Processing image {i+1}/{len(images)}: {image.id}")
                    
                    # Load image from storage
                    image_array = storage.load_image(image.image_path)
                    
                    # Validate image
                    if not cellpose_service.validate_image(image_array):
                        error_msg = f"Invalid image for segmentation: {image.id}"
                        logger.warning(error_msg)
                        results["errors"].append(error_msg)
                        results["failed_segmentations"] += 1
                        continue
                    
                    # Run segmentation
                    masks, flows = cellpose_service.segment(
                        image=image_array,
                        diameter=diameter,
                        channels=channels
                    )
                    
                    # Count cells
                    cell_count = cellpose_service.count_cells(masks)
                    results["total_cells_found"] += cell_count
                    
                    # Save segmentation mask
                    mask_path = storage.save_mask(
                        mask_data=masks,
                        segmentation_id=str(image.id)
                    )
                    
                    # Create segmentation record
                    segmentation = CellSegmentation(
                        image_id=image.id,
                        model_name=model_name,
                        model_version="3.0",
                        cell_count=cell_count,
                        mask_path=mask_path,
                        parameters={
                            "diameter": diameter,
                            "channels": channels,
                            "model_name": model_name,
                            "batch_task_id": self.request.id
                        },
                        confidence_score=0.9  # Placeholder
                    )
                    
                    db.add(segmentation)
                    db.commit()
                    db.refresh(segmentation)
                    
                    results["successful_segmentations"] += 1
                    logger.info(f"Segmented image {image.id}, found {cell_count} cells")
                    
                    # Extract features if requested
                    if extract_features and cell_count > 0:
                        try:
                            # Extract morphological features
                            morphology_features = feature_extractor.extract_morphology_features(masks)
                            
                            # Save features to database
                            features_saved = 0
                            for j, feature in enumerate(morphology_features):
                                cell_feature = CellFeature(
                                    segmentation_id=segmentation.id,
                                    cell_id=j + 1,
                                    area=feature.area,
                                    perimeter=feature.perimeter,
                                    circularity=feature.circularity,
                                    eccentricity=feature.eccentricity,
                                    solidity=feature.solidity,
                                    centroid_x=feature.centroid_x,
                                    centroid_y=feature.centroid_y,
                                    intensity_features={},
                                    texture_features={},
                                    custom_features=feature.to_dict()
                                )
                                db.add(cell_feature)
                                features_saved += 1
                            
                            db.commit()
                            results["features_extracted"] += features_saved
                            logger.info(f"Extracted features for {features_saved} cells")
                            
                        except Exception as e:
                            error_msg = f"Failed to extract features for image {image.id}: {str(e)}"
                            logger.warning(error_msg)
                            results["errors"].append(error_msg)
                    
                    results["processed_images"] += 1
                    
                    # Update progress (optional - could be used for progress tracking)
                    progress = int((i + 1) / len(images) * 100)
                    self.update_state(
                        state='PROGRESS',
                        meta={
                            'current': i + 1,
                            'total': len(images),
                            'progress': progress,
                            'status': f'Processed {i + 1}/{len(images)} images'
                        }
                    )
                    
                except Exception as e:
                    error_msg = f"Failed to process image {image.id}: {str(e)}"
                    logger.error(error_msg)
                    results["errors"].append(error_msg)
                    results["failed_segmentations"] += 1
                    continue
        
        # Calculate total processing time
        end_time = datetime.now()
        processing_time = (end_time - start_time).total_seconds()
        results["processing_time_seconds"] = processing_time
        results["completed_at"] = end_time.isoformat()
        
        logger.info(
            f"Batch segmentation completed: {results['successful_segmentations']}/{len(image_ids)} "
            f"images processed successfully in {processing_time:.1f}s"
        )
        
        return results
        
    except Exception as exc:
        error_msg = f"Batch segmentation task failed: {str(exc)}"
        logger.error(error_msg, exc_info=True)
        
        # Update results with failure info
        results["errors"].append(error_msg)
        results["failed_segmentations"] = len(image_ids)
        results["completed_at"] = datetime.now(timezone.utc).isoformat()
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            logger.error(f"Max retries reached for batch segmentation task")
            return results
        
        retry_countdown = 60 * (2 ** self.request.retries)
        logger.info(f"Retrying batch segmentation in {retry_countdown} seconds")
        self.retry(exc=exc, countdown=retry_countdown)


@celery_app.task(bind=True, max_retries=2, default_retry_delay=300, queue='low')
def cleanup_old_segmentations(self, days_old: int = 30) -> Dict[str, Any]:
    """
    Clean up old segmentation files and data.
    
    Args:
        days_old: Remove segmentations older than this many days
    
    Returns:
        Dictionary with cleanup results
    """
    try:
        from amprenta_rag.database.session import db_session
        from amprenta_rag.imaging.models import CellSegmentation
        from amprenta_rag.imaging.storage import ImageStorage
        from datetime import timedelta
        
        logger.info(f"Starting cleanup of segmentations older than {days_old} days")
        
        cutoff_date = datetime.now(timezone.utc) - timedelta(days=days_old)
        storage = ImageStorage.create_local("data/imaging")
        
        results = {
            "task_id": self.request.id,
            "cutoff_date": cutoff_date.isoformat(),
            "segmentations_found": 0,
            "segmentations_deleted": 0,
            "files_deleted": 0,
            "errors": [],
            "started_at": datetime.now(timezone.utc).isoformat()
        }
        
        with db_session() as db:
            # Find old segmentations
            old_segmentations = db.query(CellSegmentation).filter(
                CellSegmentation.created_at < cutoff_date
            ).all()
            
            results["segmentations_found"] = len(old_segmentations)
            
            for segmentation in old_segmentations:
                try:
                    # Delete mask file from storage
                    if segmentation.mask_path:
                        if storage.delete_mask(segmentation.mask_path):
                            results["files_deleted"] += 1
                        else:
                            results["errors"].append(f"Failed to delete mask file: {segmentation.mask_path}")
                    
                    # Delete database record (features will be cascade deleted)
                    db.delete(segmentation)
                    results["segmentations_deleted"] += 1
                    
                except Exception as e:
                    error_msg = f"Failed to delete segmentation {segmentation.id}: {str(e)}"
                    logger.error(error_msg)
                    results["errors"].append(error_msg)
            
            # Commit deletions
            db.commit()
        
        results["completed_at"] = datetime.now(timezone.utc).isoformat()
        
        logger.info(
            f"Cleanup completed: {results['segmentations_deleted']} segmentations "
            f"and {results['files_deleted']} files deleted"
        )
        
        return results
        
    except Exception as exc:
        error_msg = f"Cleanup task failed: {str(exc)}"
        logger.error(error_msg, exc_info=True)
        
        # Retry on failure
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": error_msg, "task_id": self.request.id}
        
        self.retry(exc=exc, countdown=300 * (2 ** self.request.retries))


@celery_app.task(queue='default')
def extract_features_for_segmentation(segmentation_id: str) -> Dict[str, Any]:
    """
    Extract features for an existing segmentation.
    
    Args:
        segmentation_id: UUID of segmentation to process
    
    Returns:
        Dictionary with extraction results
    """
    try:
        from amprenta_rag.database.session import db_session
        from amprenta_rag.imaging.models import CellSegmentation, CellFeature
        from amprenta_rag.imaging.feature_extraction import FeatureExtractor
        from amprenta_rag.imaging.storage import ImageStorage
        
        segmentation_uuid = UUID(segmentation_id)
        
        logger.info(f"Extracting features for segmentation {segmentation_id}")
        
        feature_extractor = FeatureExtractor()
        storage = ImageStorage.create_local("data/imaging")
        
        with db_session() as db:
            # Get segmentation
            segmentation = db.query(CellSegmentation).filter(
                CellSegmentation.id == segmentation_uuid
            ).first()
            
            if not segmentation:
                return {"status": "failed", "error": "Segmentation not found"}
            
            # Load mask
            masks = storage.load_mask(segmentation.mask_path)
            
            # Extract features
            morphology_features = feature_extractor.extract_morphology_features(masks)
            
            # Save features
            features_saved = 0
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
                    intensity_features={},
                    texture_features={},
                    custom_features=feature.to_dict()
                )
                db.add(cell_feature)
                features_saved += 1
            
            db.commit()
            
            logger.info(f"Extracted features for {features_saved} cells")
            
            return {
                "status": "success",
                "segmentation_id": segmentation_id,
                "features_extracted": features_saved,
                "cell_count": segmentation.cell_count
            }
            
    except Exception as e:
        error_msg = f"Feature extraction failed: {str(e)}"
        logger.error(error_msg, exc_info=True)
        return {"status": "failed", "error": error_msg}
