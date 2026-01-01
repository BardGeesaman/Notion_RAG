"""Unit tests for imaging API endpoints."""

from __future__ import annotations

import io
import uuid
from unittest.mock import MagicMock, patch, Mock
from datetime import datetime

import pytest
from fastapi.testclient import TestClient
from PIL import Image

from amprenta_rag.api.main import app
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.models.chemistry import HTSWell, HTSPlate, HTSCampaign


class TestImagingAPI:
    """Test imaging API endpoints."""

    def setup_method(self):
        """Set up test client and mock dependencies."""
        self.client = TestClient(app)
        
        # Mock authentication
        self.mock_user = MagicMock()
        self.mock_user.id = uuid.uuid4()
        
        # Override dependencies
        app.dependency_overrides[get_current_user] = lambda: self.mock_user
        app.dependency_overrides[get_db] = self.mock_db_session
        app.dependency_overrides[get_cellpose_service] = self.mock_cellpose_service
        app.dependency_overrides[get_feature_extractor] = self.mock_feature_extractor
        app.dependency_overrides[get_image_storage] = self.mock_image_storage

    def teardown_method(self):
        """Clean up after each test."""
        app.dependency_overrides.clear()

    def mock_db_session(self):
        """Mock database session."""
        session = MagicMock()
        return session

    def mock_cellpose_service(self):
        """Mock CellPose service."""
        service = MagicMock()
        service.validate_image.return_value = True
        service.segment.return_value = (
            [[0, 1, 1], [0, 1, 1], [2, 2, 0]],  # masks
            [[[0, 0], [1, 1], [1, 1]], [[0, 0], [1, 1], [1, 1]], [[2, 2], [2, 2], [0, 0]]]  # flows
        )
        service.count_cells.return_value = 2
        return service

    def mock_feature_extractor(self):
        """Mock feature extractor."""
        extractor = MagicMock()
        
        # Mock morphology features
        mock_feature = MagicMock()
        mock_feature.area = 100.0
        mock_feature.perimeter = 40.0
        mock_feature.circularity = 0.8
        mock_feature.eccentricity = 0.5
        mock_feature.solidity = 0.9
        mock_feature.centroid_x = 50.0
        mock_feature.centroid_y = 60.0
        mock_feature.to_dict.return_value = {
            "area": 100.0, "perimeter": 40.0, "circularity": 0.8
        }
        
        extractor.extract_morphology_features.return_value = [mock_feature, mock_feature]
        
        # Mock aggregated features
        mock_aggregated = MagicMock()
        mock_aggregated.to_dict.return_value = {
            "cell_count": 2,
            "morphology_stats": {"area": {"mean": 100.0, "std": 10.0}},
            "intensity_stats": {}
        }
        extractor.aggregate_to_well.return_value = mock_aggregated
        extractor.calculate_well_metrics.return_value = {"cell_count": 2, "mean_cell_area": 100.0}
        
        return extractor

    def mock_image_storage(self):
        """Mock image storage."""
        storage = MagicMock()
        storage.save_image.return_value = "images/test/DAPI_z000_t000.png"
        storage.save_mask.return_value = "masks/test_mask.npy"
        storage.load_image.return_value = [[100, 150, 200], [120, 160, 180], [110, 140, 190]]
        storage.load_mask.return_value = [[0, 1, 1], [0, 1, 1], [2, 2, 0]]
        return storage

    def create_test_image_file(self):
        """Create a test image file for upload."""
        # Create a simple test image
        img = Image.new('L', (100, 100), color=128)
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)
        return img_bytes

    def test_upload_image(self):
        """Test uploading an image without well association."""
        # Mock database session
        mock_db = MagicMock()
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Create test image file
        test_image = self.create_test_image_file()
        
        # Upload image
        response = self.client.post(
            "/api/v1/imaging/upload",
            files={"file": ("test.png", test_image, "image/png")},
            data={
                "channel": "DAPI",
                "z_slice": 0,
                "timepoint": 0,
                "pixel_size_um": 0.325
            }
        )
        
        assert response.status_code == 201
        data = response.json()
        assert "image_id" in data
        assert data["channel"] == "DAPI"
        assert data["width"] == 100
        assert data["height"] == 100
        assert data["well_id"] is None
        assert "message" in data

    def test_upload_image_with_well(self):
        """Test uploading an image with well association."""
        # Mock database session and well lookup
        mock_db = MagicMock()
        mock_well = MagicMock()
        mock_well.id = uuid.uuid4()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_well
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Create test image file
        test_image = self.create_test_image_file()
        
        well_id = str(mock_well.id)
        
        # Upload image with well association
        response = self.client.post(
            "/api/v1/imaging/upload",
            files={"file": ("test.png", test_image, "image/png")},
            data={
                "well_id": well_id,
                "channel": "GFP",
                "z_slice": 1,
                "timepoint": 2
            }
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["well_id"] == well_id
        assert data["channel"] == "GFP"

    def test_segment_image(self):
        """Test segmenting an image."""
        # Mock database session and image lookup
        mock_db = MagicMock()
        mock_image = MagicMock()
        mock_image.id = uuid.uuid4()
        mock_image.image_path = "images/test.png"
        mock_db.query.return_value.filter.return_value.first.return_value = mock_image
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Segment image
        response = self.client.post(
            "/api/v1/imaging/segment",
            json={
                "image_id": str(mock_image.id),
                "model_name": "cyto",
                "diameter": 30.0,
                "channels": [0, 0],
                "extract_features": True
            }
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "segmentation_id" in data
        assert data["image_id"] == str(mock_image.id)
        assert data["cell_count"] == 2
        assert data["model_name"] == "cyto"
        assert data["features_extracted"] is True

    def test_segment_invalid_image(self):
        """Test segmenting a non-existent image returns 404."""
        # Mock database session with no image found
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Try to segment non-existent image
        response = self.client.post(
            "/api/v1/imaging/segment",
            json={
                "image_id": str(uuid.uuid4()),
                "model_name": "cyto"
            }
        )
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    @patch('amprenta_rag.jobs.tasks.imaging.process_batch_segmentation')
    def test_segment_batch_queues_task(self, mock_task):
        """Test that batch segmentation queues a Celery task."""
        # Mock Celery task
        mock_task.delay.return_value.id = "test-task-123"
        
        # Mock database session with existing images
        mock_db = MagicMock()
        mock_image1 = MagicMock()
        mock_image1.id = uuid.uuid4()
        mock_image2 = MagicMock()
        mock_image2.id = uuid.uuid4()
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_image1, mock_image2]
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Queue batch segmentation
        response = self.client.post(
            "/api/v1/imaging/segment-batch",
            json={
                "image_ids": [str(mock_image1.id), str(mock_image2.id)],
                "model_name": "cyto",
                "extract_features": True
            }
        )
        
        assert response.status_code == 202
        data = response.json()
        assert data["task_id"] == "test-task-123"
        assert data["image_count"] == 2
        assert data["status"] == "queued"
        
        # Verify task was called
        mock_task.delay.assert_called_once()

    def test_get_image_metadata(self):
        """Test getting image metadata."""
        # Mock database session and image
        mock_db = MagicMock()
        mock_image = MagicMock()
        mock_image.id = uuid.uuid4()
        mock_image.well_id = uuid.uuid4()
        mock_image.channel = "DAPI"
        mock_image.z_slice = 0
        mock_image.timepoint = 0
        mock_image.width = 1024
        mock_image.height = 1024
        mock_image.bit_depth = 16
        mock_image.pixel_size_um = 0.325
        mock_image.image_path = "images/test.tiff"
        mock_image.image_metadata = {"exposure_ms": 100}
        mock_image.acquired_at = datetime.now()
        mock_image.created_at = datetime.now()
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_image
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Get image metadata
        response = self.client.get(f"/api/v1/imaging/images/{mock_image.id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["image_id"] == str(mock_image.id)
        assert data["channel"] == "DAPI"
        assert data["width"] == 1024
        assert data["height"] == 1024

    def test_get_segmentation(self):
        """Test getting segmentation results."""
        # Mock database session and segmentation
        mock_db = MagicMock()
        mock_segmentation = MagicMock()
        mock_segmentation.id = uuid.uuid4()
        mock_segmentation.image_id = uuid.uuid4()
        mock_segmentation.model_name = "cyto"
        mock_segmentation.model_version = "3.0"
        mock_segmentation.cell_count = 150
        mock_segmentation.mask_path = "masks/test_mask.npy"
        mock_segmentation.parameters = {"diameter": 30}
        mock_segmentation.confidence_score = 0.95
        mock_segmentation.created_at = datetime.now()
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_segmentation
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Get segmentation
        response = self.client.get(f"/api/v1/imaging/segmentations/{mock_segmentation.id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["segmentation_id"] == str(mock_segmentation.id)
        assert data["cell_count"] == 150
        assert data["model_name"] == "cyto"

    def test_get_features(self):
        """Test getting cell features."""
        # Mock database session
        mock_db = MagicMock()
        
        # Mock segmentation
        mock_segmentation = MagicMock()
        mock_segmentation.id = uuid.uuid4()
        
        # Mock features
        mock_feature1 = MagicMock()
        mock_feature1.cell_id = 1
        mock_feature1.area = 100.0
        mock_feature1.perimeter = 40.0
        mock_feature1.circularity = 0.8
        mock_feature1.eccentricity = 0.5
        mock_feature1.solidity = 0.9
        mock_feature1.centroid_x = 50.0
        mock_feature1.centroid_y = 60.0
        mock_feature1.intensity_features = {}
        mock_feature1.texture_features = {}
        mock_feature1.custom_features = {}
        
        mock_feature2 = MagicMock()
        mock_feature2.cell_id = 2
        mock_feature2.area = 120.0
        mock_feature2.perimeter = 45.0
        mock_feature2.circularity = 0.75
        mock_feature2.eccentricity = 0.6
        mock_feature2.solidity = 0.85
        mock_feature2.centroid_x = 80.0
        mock_feature2.centroid_y = 90.0
        mock_feature2.intensity_features = {}
        mock_feature2.texture_features = {}
        mock_feature2.custom_features = {}
        
        # Set up query mocks
        mock_db.query.return_value.filter.return_value.first.return_value = mock_segmentation
        mock_db.query.return_value.filter.return_value.all.return_value = [mock_feature1, mock_feature2]
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Get features
        response = self.client.get(f"/api/v1/imaging/features/{mock_segmentation.id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["segmentation_id"] == str(mock_segmentation.id)
        assert data["cell_count"] == 2
        assert len(data["features"]) == 2
        assert data["features"][0]["cell_id"] == 1
        assert data["features"][0]["area"] == 100.0

    def test_get_well_summary(self):
        """Test getting well summary with aggregated features."""
        # Mock database session
        mock_db = MagicMock()
        
        # Mock well
        mock_well = MagicMock()
        mock_well.id = uuid.uuid4()
        
        # Mock images
        mock_image1 = MagicMock()
        mock_image1.id = uuid.uuid4()
        mock_image1.channel = "DAPI"
        mock_image2 = MagicMock()
        mock_image2.id = uuid.uuid4()
        mock_image2.channel = "GFP"
        
        # Mock segmentations
        mock_seg1 = MagicMock()
        mock_seg1.id = uuid.uuid4()
        mock_seg1.cell_count = 50
        mock_seg2 = MagicMock()
        mock_seg2.id = uuid.uuid4()
        mock_seg2.cell_count = 60
        
        # Mock features
        mock_features = [MagicMock() for _ in range(10)]
        for i, feature in enumerate(mock_features):
            feature.area = 100.0 + i * 10
            feature.perimeter = 40.0 + i * 2
            feature.circularity = 0.8
            feature.eccentricity = 0.5
            feature.solidity = 0.9
            feature.centroid_x = 50.0
            feature.centroid_y = 60.0
        
        # Set up query chain
        queries = [
            mock_well,  # Well query
            [mock_image1, mock_image2],  # Images query
            [mock_seg1, mock_seg2],  # Segmentations query
            mock_features  # Features query
        ]
        
        def side_effect(*args, **kwargs):
            return Mock(filter=Mock(return_value=Mock(
                first=Mock(return_value=queries[0] if len(queries) == 1 else None),
                all=Mock(return_value=queries.pop(0) if queries else [])
            )))
        
        mock_db.query.side_effect = side_effect
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Get well summary
        response = self.client.get(f"/api/v1/imaging/wells/{mock_well.id}/summary")
        
        assert response.status_code == 200
        data = response.json()
        assert data["well_id"] == str(mock_well.id)
        assert data["image_count"] == 2
        assert data["total_cell_count"] == 110  # 50 + 60
        assert "DAPI" in data["channels"]
        assert "GFP" in data["channels"]

    def test_get_well_summary_no_images(self):
        """Test getting well summary for empty well."""
        # Mock database session
        mock_db = MagicMock()
        
        # Mock well exists but no images
        mock_well = MagicMock()
        mock_well.id = uuid.uuid4()
        
        # Set up queries to return well but no images
        mock_db.query.return_value.filter.return_value.first.return_value = mock_well
        mock_db.query.return_value.filter.return_value.all.return_value = []
        app.dependency_overrides[get_db] = lambda: mock_db
        
        # Get well summary
        response = self.client.get(f"/api/v1/imaging/wells/{mock_well.id}/summary")
        
        assert response.status_code == 200
        data = response.json()
        assert data["well_id"] == str(mock_well.id)
        assert data["image_count"] == 0
        assert data["total_cell_count"] == 0
        assert data["channels"] == []
        assert data["aggregated_features"] == {}


# Import required dependencies for testing
try:
    from amprenta_rag.api.dependencies import get_current_user, get_db
    from amprenta_rag.api.routers.imaging import (
        get_cellpose_service, get_feature_extractor, get_image_storage
    )
except ImportError:
    # Handle missing imports for test environment
    def get_current_user():
        pass
    
    def get_db():
        pass
    
    def get_cellpose_service():
        pass
    
    def get_feature_extractor():
        pass
    
    def get_image_storage():
        pass
