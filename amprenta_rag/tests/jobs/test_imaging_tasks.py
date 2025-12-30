"""Tests for imaging Celery tasks."""

import os
from datetime import datetime, timezone
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.jobs.tasks.imaging import (
    process_batch_segmentation,
    cleanup_old_segmentations,
    extract_features_for_segmentation
)


class TestImagingTasks:
    """Test imaging Celery tasks."""
    
    def setup_method(self):
        os.environ["CELERY_TASK_ALWAYS_EAGER"] = "true"
    
    def teardown_method(self):
        os.environ.pop("CELERY_TASK_ALWAYS_EAGER", None)


class TestBatchSegmentation(TestImagingTasks):
    """Test batch segmentation task."""
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.models.CellFeature')
    @patch('amprenta_rag.imaging.models.CellSegmentation')
    @patch('amprenta_rag.imaging.models.MicroscopyImage')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.feature_extraction.FeatureExtractor')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_success(self, mock_storage_class, mock_extractor_class, 
                                       mock_cellpose_class, mock_image_class,
                                       mock_segmentation_class, mock_feature_class, mock_db_session):
        """Test successful batch segmentation with multiple images."""
        image_id1 = uuid4()
        image_id2 = uuid4()
        image_ids = [str(image_id1), str(image_id2)]
        
        # Mock images
        mock_image1 = MagicMock()
        mock_image1.id = image_id1
        mock_image1.image_path = "/path/to/image1.tif"
        mock_image2 = MagicMock()
        mock_image2.id = image_id2
        mock_image2.image_path = "/path/to/image2.tif"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image1, mock_image2]
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.return_value = MagicMock()  # Mock image array
        mock_storage.save_mask.return_value = "/path/to/mask.npy"
        
        # Mock CellPose service
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = True
        mock_cellpose.segment.return_value = (MagicMock(), MagicMock())  # masks, flows
        mock_cellpose.count_cells.return_value = 150
        
        # Mock feature extractor
        mock_extractor = MagicMock()
        mock_extractor_class.return_value = mock_extractor
        mock_feature = MagicMock()
        mock_feature.area = 100.0
        mock_feature.perimeter = 35.0
        mock_feature.circularity = 0.8
        mock_feature.eccentricity = 0.3
        mock_feature.solidity = 0.9
        mock_feature.centroid_x = 50.0
        mock_feature.centroid_y = 75.0
        mock_feature.to_dict.return_value = {"custom_metric": 1.5}
        mock_extractor.extract_morphology_features.return_value = [mock_feature] * 150
        
        # Execute task
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(
                image_ids=image_ids,
                model_name="cyto",
                diameter=30.0,
                channels=[0, 0],
                extract_features=True,
                user_id="test_user"
            )
        
        # Debug: Print result if test fails
        if result["failed_segmentations"] != 0:
            print(f"DEBUG: Result = {result}")
        
        # Verify result
        assert result["total_images"] == 2
        assert result["processed_images"] == 2
        assert result["successful_segmentations"] == 2
        assert result["failed_segmentations"] == 0
        assert result["total_cells_found"] == 300  # 150 per image
        assert result["features_extracted"] == 300  # 150 per image
        assert result["user_id"] == "test_user"
        assert len(result["errors"]) == 0
        
        # Verify services were called correctly
        mock_cellpose_class.assert_called_once_with(
            model_type="cyto",
            gpu=True,
            tile_size=1024,
            overlap=128
        )
        assert mock_cellpose.segment.call_count == 2
        assert mock_extractor.extract_morphology_features.call_count == 2
    
    @patch('amprenta_rag.database.session.db_session')
    def test_batch_segmentation_image_not_found(self, mock_db_session):
        """Test batch segmentation when images are missing."""
        image_id1 = uuid4()
        image_id2 = uuid4()
        image_ids = [str(image_id1), str(image_id2)]
        
        # Mock database session with only one image found
        mock_image1 = MagicMock()
        mock_image1.id = image_id1
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image1]  # Missing one image
        
        # Execute task
        result = process_batch_segmentation(image_ids=image_ids)
        
        # Verify result shows error for missing image
        assert result["total_images"] == 2
        assert result["processed_images"] == 0
        assert result["successful_segmentations"] == 0
        assert result["failed_segmentations"] == 0
        assert len(result["errors"]) == 1
        assert "Images not found" in result["errors"][0]
        assert str(image_id2) in result["errors"][0]
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_invalid_image(self, mock_storage_class, mock_cellpose_class, mock_db_session):
        """Test batch segmentation when image validation fails."""
        image_id = uuid4()
        image_ids = [str(image_id)]
        
        # Mock image
        mock_image = MagicMock()
        mock_image.id = image_id
        mock_image.image_path = "/path/to/invalid.tif"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image]
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.return_value = MagicMock()  # Mock image array
        
        # Mock CellPose service with validation failure
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = False  # Invalid image
        
        # Execute task
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(image_ids=image_ids)
        
        # Verify result shows validation error
        assert result["total_images"] == 1
        assert result["processed_images"] == 1
        assert result["successful_segmentations"] == 0
        assert result["failed_segmentations"] == 1
        assert len(result["errors"]) == 1
        assert "Invalid image for segmentation" in result["errors"][0]
        
        # Verify segmentation was not called
        mock_cellpose.segment.assert_not_called()
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.feature_extraction.FeatureExtractor')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_partial_success(self, mock_storage_class, mock_extractor_class,
                                               mock_cellpose_class, mock_db_session):
        """Test batch segmentation when some images fail but processing continues."""
        image_id1 = uuid4()
        image_id2 = uuid4()
        image_ids = [str(image_id1), str(image_id2)]
        
        # Mock images
        mock_image1 = MagicMock()
        mock_image1.id = image_id1
        mock_image1.image_path = "/path/to/image1.tif"
        mock_image2 = MagicMock()
        mock_image2.id = image_id2
        mock_image2.image_path = "/path/to/image2.tif"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image1, mock_image2]
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.side_effect = [MagicMock(), Exception("File corrupted")]  # Second image fails
        mock_storage.save_mask.return_value = "/path/to/mask.npy"
        
        # Mock CellPose service
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = True
        mock_cellpose.segment.return_value = (MagicMock(), MagicMock())
        mock_cellpose.count_cells.return_value = 100
        
        # Mock feature extractor
        mock_extractor = MagicMock()
        mock_extractor_class.return_value = mock_extractor
        mock_feature = MagicMock()
        mock_feature.area = 80.0
        mock_feature.to_dict.return_value = {}
        mock_extractor.extract_morphology_features.return_value = [mock_feature] * 100
        
        # Execute task
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(image_ids=image_ids, extract_features=True)
        
        # Verify partial success
        assert result["total_images"] == 2
        assert result["processed_images"] == 2
        assert result["successful_segmentations"] == 1  # Only first image succeeded
        assert result["failed_segmentations"] == 1     # Second image failed
        assert result["total_cells_found"] == 100      # Only from first image
        assert result["features_extracted"] == 100     # Only from first image
        assert len(result["errors"]) == 1
        assert "File corrupted" in result["errors"][0]
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.feature_extraction.FeatureExtractor')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_feature_extraction(self, mock_storage_class, mock_extractor_class,
                                                   mock_cellpose_class, mock_db_session):
        """Test batch segmentation validates feature extraction."""
        image_id = uuid4()
        image_ids = [str(image_id)]
        
        # Mock image
        mock_image = MagicMock()
        mock_image.id = image_id
        mock_image.image_path = "/path/to/image.tif"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image]
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.return_value = MagicMock()
        mock_storage.save_mask.return_value = "/path/to/mask.npy"
        
        # Mock CellPose service
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = True
        mock_cellpose.segment.return_value = (MagicMock(), MagicMock())
        mock_cellpose.count_cells.return_value = 50
        
        # Mock feature extractor with specific features
        mock_extractor = MagicMock()
        mock_extractor_class.return_value = mock_extractor
        mock_features = []
        for i in range(50):
            mock_feature = MagicMock()
            mock_feature.area = 100.0 + i
            mock_feature.perimeter = 35.0 + i * 0.5
            mock_feature.circularity = 0.8
            mock_feature.eccentricity = 0.3
            mock_feature.solidity = 0.9
            mock_feature.centroid_x = 50.0 + i
            mock_feature.centroid_y = 75.0 + i
            mock_feature.to_dict.return_value = {"custom": i}
            mock_features.append(mock_feature)
        mock_extractor.extract_morphology_features.return_value = mock_features
        
        # Execute task with feature extraction
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(image_ids=image_ids, extract_features=True)
        
        # Verify feature extraction was performed
        assert result["features_extracted"] == 50
        mock_extractor.extract_morphology_features.assert_called_once()
        
        # Verify CellFeature records were added to database
        assert mock_db.add.call_count >= 50  # At least 50 feature records + segmentation
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_retry_on_failure(self, mock_storage_class, mock_cellpose_class, mock_db_session):
        """Test batch segmentation handles failure and re-raises for retry."""
        image_ids = [str(uuid4())]
        
        # Mock database session to raise exception
        mock_db_session.side_effect = Exception("Database connection failed")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            process_batch_segmentation(image_ids=image_ids)
        
        # Verify correct error
        assert "Database connection failed" in str(exc_info.value)
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_progress_tracking(self, mock_storage_class, mock_cellpose_class, mock_db_session):
        """Test batch segmentation results dict has correct counts."""
        image_id1 = uuid4()
        image_id2 = uuid4()
        image_id3 = uuid4()
        image_ids = [str(image_id1), str(image_id2), str(image_id3)]
        
        # Mock images
        mock_images = []
        for i, img_id in enumerate([image_id1, image_id2, image_id3]):
            mock_image = MagicMock()
            mock_image.id = img_id
            mock_image.image_path = f"/path/to/image{i}.tif"
            mock_images.append(mock_image)
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = mock_images
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.return_value = MagicMock()
        mock_storage.save_mask.return_value = "/path/to/mask.npy"
        
        # Mock CellPose service
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = True
        mock_cellpose.segment.return_value = (MagicMock(), MagicMock())
        mock_cellpose.count_cells.side_effect = [75, 100, 125]  # Different cell counts
        
        # Execute task
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(image_ids=image_ids, extract_features=False)
        
        # Verify progress tracking
        assert result["total_images"] == 3
        assert result["processed_images"] == 3
        assert result["successful_segmentations"] == 3
        assert result["failed_segmentations"] == 0
        assert result["total_cells_found"] == 300  # 75 + 100 + 125
        assert result["features_extracted"] == 0   # Features not extracted
        assert "started_at" in result
        assert "completed_at" in result
        assert "processing_time_seconds" in result
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.cellpose_service.CellPoseService')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_batch_segmentation_db_records(self, mock_storage_class, mock_cellpose_class, mock_db_session):
        """Test batch segmentation creates CellSegmentation records."""
        image_id = uuid4()
        image_ids = [str(image_id)]
        
        # Mock image
        mock_image = MagicMock()
        mock_image.id = image_id
        mock_image.image_path = "/path/to/image.tif"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image]
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_image.return_value = MagicMock()
        mock_storage.save_mask.return_value = "/path/to/mask.npy"
        
        # Mock CellPose service
        mock_cellpose = MagicMock()
        mock_cellpose_class.return_value = mock_cellpose
        mock_cellpose.validate_image.return_value = True
        mock_cellpose.segment.return_value = (MagicMock(), MagicMock())
        mock_cellpose.count_cells.return_value = 80
        
        # Execute task
        with patch.object(process_batch_segmentation, 'update_state'):
            result = process_batch_segmentation(image_ids=image_ids, model_name="nuclei", diameter=25.0)
        
        # Verify segmentation record creation
        assert result["successful_segmentations"] == 1
        assert result["total_cells_found"] == 80
        
        # Verify CellSegmentation was added to database
        mock_db.add.assert_called()
        mock_db.commit.assert_called()
        mock_db.refresh.assert_called()


class TestCleanupSegmentations(TestImagingTasks):
    """Test cleanup segmentations task."""
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_cleanup_success(self, mock_storage_class, mock_db_session):
        """Test successful cleanup of old segmentations."""
        # Mock old segmentations
        mock_seg1 = MagicMock()
        mock_seg1.id = uuid4()
        mock_seg1.mask_path = "/path/to/old_mask1.npy"
        mock_seg2 = MagicMock()
        mock_seg2.id = uuid4()
        mock_seg2.mask_path = "/path/to/old_mask2.npy"
        old_segmentations = [mock_seg1, mock_seg2]
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = old_segmentations
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.delete_mask.return_value = True
        
        # Execute task
        result = cleanup_old_segmentations(days_old=30)
        
        # Verify result
        assert result["segmentations_found"] == 2
        assert result["segmentations_deleted"] == 2
        assert result["files_deleted"] == 2
        assert len(result["errors"]) == 0
        
        # Verify storage cleanup was called
        mock_storage.delete_mask.assert_any_call("/path/to/old_mask1.npy")
        mock_storage.delete_mask.assert_any_call("/path/to/old_mask2.npy")
        
        # Verify database deletions
        mock_db.delete.assert_any_call(mock_seg1)
        mock_db.delete.assert_any_call(mock_seg2)
        mock_db.commit.assert_called_once()
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_cleanup_no_old_data(self, mock_storage_class, mock_db_session):
        """Test cleanup when no old segmentations exist."""
        # Mock database session with empty result
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = []
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        
        # Execute task
        result = cleanup_old_segmentations(days_old=60)
        
        # Verify result
        assert result["segmentations_found"] == 0
        assert result["segmentations_deleted"] == 0
        assert result["files_deleted"] == 0
        assert len(result["errors"]) == 0
        
        # Verify no deletion calls were made
        mock_storage.delete_mask.assert_not_called()
        mock_db.delete.assert_not_called()
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_cleanup_file_deletion(self, mock_storage_class, mock_db_session):
        """Test cleanup verifies file deletion."""
        # Mock old segmentation
        mock_segmentation = MagicMock()
        mock_segmentation.id = uuid4()
        mock_segmentation.mask_path = "/path/to/mask.npy"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_segmentation]
        
        # Mock storage with file deletion failure
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.delete_mask.return_value = False  # File deletion failed
        
        # Execute task
        result = cleanup_old_segmentations(days_old=15)
        
        # Verify result shows file deletion error
        assert result["segmentations_found"] == 1
        assert result["segmentations_deleted"] == 1  # DB record still deleted
        assert result["files_deleted"] == 0         # File deletion failed
        assert len(result["errors"]) == 1
        assert "Failed to delete mask file" in result["errors"][0]
    
    @patch('amprenta_rag.database.session.db_session')
    def test_cleanup_retry_on_failure(self, mock_db_session):
        """Test cleanup task handles failure and re-raises for retry."""
        # Mock database session to raise exception
        mock_db_session.side_effect = Exception("Database error")
        
        # Execute task - should raise exception for retry
        with pytest.raises(Exception) as exc_info:
            cleanup_old_segmentations(days_old=30)
        
        # Verify correct error
        assert "Database error" in str(exc_info.value)


class TestExtractFeatures(TestImagingTasks):
    """Test extract features task."""
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.models.CellFeature')
    @patch('amprenta_rag.imaging.models.CellSegmentation')
    @patch('amprenta_rag.imaging.feature_extraction.FeatureExtractor')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_extract_features_success(self, mock_storage_class, mock_extractor_class, 
                                     mock_segmentation_class, mock_feature_class, mock_db_session):
        """Test successful feature extraction."""
        segmentation_id = uuid4()
        
        # Mock segmentation
        mock_segmentation = MagicMock()
        mock_segmentation.id = segmentation_id
        mock_segmentation.mask_path = "/path/to/mask.npy"
        mock_segmentation.cell_count = 25
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = mock_segmentation
        
        # Mock storage
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_mask.return_value = MagicMock()  # Mock mask data
        
        # Mock feature extractor
        mock_extractor = MagicMock()
        mock_extractor_class.return_value = mock_extractor
        mock_features = []
        for i in range(25):
            mock_feature = MagicMock()
            mock_feature.area = 100.0 + i
            mock_feature.perimeter = 35.0
            mock_feature.circularity = 0.8
            mock_feature.eccentricity = 0.3
            mock_feature.solidity = 0.9
            mock_feature.centroid_x = 50.0
            mock_feature.centroid_y = 75.0
            mock_feature.to_dict.return_value = {"custom": i}
            mock_features.append(mock_feature)
        mock_extractor.extract_morphology_features.return_value = mock_features
        
        # Execute task
        result = extract_features_for_segmentation(str(segmentation_id))
        
        # Verify result
        assert result["status"] == "success"
        assert result["segmentation_id"] == str(segmentation_id)
        assert result["features_extracted"] == 25
        assert result["cell_count"] == 25
        
        # Verify feature extraction was performed
        mock_extractor.extract_morphology_features.assert_called_once()
        mock_storage.load_mask.assert_called_once_with("/path/to/mask.npy")
        
        # Verify CellFeature records were added
        assert mock_db.add.call_count == 25
        mock_db.commit.assert_called_once()
    
    @patch('amprenta_rag.database.session.db_session')
    def test_extract_features_not_found(self, mock_db_session):
        """Test feature extraction when segmentation doesn't exist."""
        segmentation_id = uuid4()
        
        # Mock database session with no segmentation found
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        # Execute task
        result = extract_features_for_segmentation(str(segmentation_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert result["error"] == "Segmentation not found"
    
    @patch('amprenta_rag.database.session.db_session')
    @patch('amprenta_rag.imaging.storage.ImageStorage')
    def test_extract_features_mask_error(self, mock_storage_class, mock_db_session):
        """Test feature extraction when mask loading fails."""
        segmentation_id = uuid4()
        
        # Mock segmentation
        mock_segmentation = MagicMock()
        mock_segmentation.id = segmentation_id
        mock_segmentation.mask_path = "/path/to/missing_mask.npy"
        
        # Mock database session
        mock_db = MagicMock()
        mock_db_session.return_value.__enter__.return_value = mock_db
        mock_db.query.return_value.filter.return_value.first.return_value = mock_segmentation
        
        # Mock storage with mask loading failure
        mock_storage = MagicMock()
        mock_storage_class.create_local.return_value = mock_storage
        mock_storage.load_mask.side_effect = FileNotFoundError("Mask file not found")
        
        # Execute task
        result = extract_features_for_segmentation(str(segmentation_id))
        
        # Verify result
        assert result["status"] == "failed"
        assert "Mask file not found" in result["error"]
