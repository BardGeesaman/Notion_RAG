"""Integration tests for imaging API with real database."""

import pytest
import io
import tempfile
from unittest.mock import patch, MagicMock
from uuid import uuid4
from datetime import datetime
from PIL import Image

from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.models.chemistry import HTSWell, HTSPlate, HTSCampaign


@pytest.mark.integration
class TestImagingAPIIntegration:
    """Integration tests for imaging API endpoints."""

    @pytest.fixture
    def test_campaign(self, db_session):
        """Create a test HTS campaign."""
        campaign = HTSCampaign(
            id=uuid4(),
            name=f"Test Campaign {uuid4().hex[:8]}",
            description="Test campaign for imaging"
        )
        db_session.add(campaign)
        db_session.commit()
        db_session.refresh(campaign)
        return campaign

    @pytest.fixture
    def test_plate(self, db_session, test_campaign):
        """Create a test HTS plate."""
        plate = HTSPlate(
            id=uuid4(),
            name=f"Test Plate {uuid4().hex[:8]}",
            campaign_id=test_campaign.id,
            plate_format=384
        )
        db_session.add(plate)
        db_session.commit()
        db_session.refresh(plate)
        return plate

    @pytest.fixture
    def test_well(self, db_session, test_plate):
        """Create a test HTS well."""
        well = HTSWell(
            id=uuid4(),
            plate_id=test_plate.id,
            row=1,
            column=1,
            well_id="A01"
        )
        db_session.add(well)
        db_session.commit()
        db_session.refresh(well)
        return well

    def create_test_image_file(self):
        """Create a test image file for upload."""
        img = Image.new('L', (100, 100), color=128)
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)
        return img_bytes

    @patch('amprenta_rag.api.routers.imaging.get_image_storage')
    @patch('amprenta_rag.api.routers.imaging.get_cellpose_service')
    def test_upload_image_integration(self, mock_cellpose, mock_storage, 
                                    integration_client, db_session, timed_request):
        """Test uploading an image without well association."""
        # Mock external services (storage and cellpose)
        mock_storage.return_value.save_image.return_value = "images/test/DAPI_z000_t000.png"
        mock_cellpose.return_value.validate_image.return_value = True
        
        # Create test image file
        test_image = self.create_test_image_file()
        
        # Upload image via API
        files = {"file": ("test.png", test_image, "image/png")}
        data = {
            "channel": "DAPI",
            "z_slice": "0",
            "timepoint": "0", 
            "pixel_size_um": "0.325"
        }
        
        response = integration_client.post("/api/v1/imaging/upload", files=files, data=data)
        
        assert response.status_code == 201
        response_data = response.json()
        assert "image_id" in response_data
        assert response_data["channel"] == "DAPI"
        assert response_data["width"] == 100
        assert response_data["height"] == 100
        assert response_data["well_id"] is None
        
        # Verify image was created in database
        image_id = response_data["image_id"]
        image = db_session.query(MicroscopyImage).filter_by(id=image_id).first()
        assert image is not None
        assert image.channel == "DAPI"
        assert image.width == 100
        assert image.height == 100
        assert image.pixel_size_um == 0.325

    @patch('amprenta_rag.api.routers.imaging.get_image_storage')
    @patch('amprenta_rag.api.routers.imaging.get_cellpose_service')
    def test_upload_image_with_well_integration(self, mock_cellpose, mock_storage,
                                              integration_client, db_session, test_well, timed_request):
        """Test uploading an image with well association."""
        # Mock external services
        mock_storage.return_value.save_image.return_value = "images/test/DAPI_z000_t000.png"
        mock_cellpose.return_value.validate_image.return_value = True
        
        # Create test image file
        test_image = self.create_test_image_file()
        
        # Upload image with well association
        files = {"file": ("test.png", test_image, "image/png")}
        data = {
            "channel": "DAPI",
            "z_slice": "0",
            "timepoint": "0",
            "pixel_size_um": "0.325",
            "well_id": str(test_well.id)
        }
        
        response = integration_client.post("/api/v1/imaging/upload", files=files, data=data)
        
        assert response.status_code == 201
        response_data = response.json()
        assert response_data["well_id"] == str(test_well.id)
        
        # Verify image was linked to well in database
        image_id = response_data["image_id"]
        image = db_session.query(MicroscopyImage).filter_by(id=image_id).first()
        assert image is not None
        assert image.well_id == test_well.id

    @patch('amprenta_rag.api.routers.imaging.get_cellpose_service')
    @patch('amprenta_rag.api.routers.imaging.get_image_storage')
    def test_segment_image_integration(self, mock_storage, mock_cellpose,
                                     integration_client, db_session, timed_request):
        """Test image segmentation with real database."""
        # Create real image in database
        image = MicroscopyImage(
            id=uuid4(),
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=100,
            height=100,
            pixel_size_um=0.325,
            image_path="images/test/DAPI_z000_t000.png",
            created_at=datetime.utcnow()
        )
        db_session.add(image)
        db_session.commit()
        db_session.refresh(image)
        
        # Mock external services
        mock_storage.return_value.load_image.return_value = [[100, 150, 200], [120, 160, 180]]
        mock_storage.return_value.save_mask.return_value = "masks/test_mask.npy"
        mock_cellpose.return_value.segment.return_value = (
            [[0, 1, 1], [0, 1, 1], [2, 2, 0]],  # masks
            [[[0, 0], [1, 1], [1, 1]], [[0, 0], [1, 1], [1, 1]], [[2, 2], [2, 2], [0, 0]]]  # flows
        )
        mock_cellpose.return_value.count_cells.return_value = 2
        
        response, benchmark = timed_request(
            "POST",
            f"/api/v1/imaging/{image.id}/segment",
            "test_segment_image",
            json={"model": "cyto", "diameter": 30}
        )
        
        assert response.status_code == 201
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert "segmentation_id" in response_data
        assert response_data["cell_count"] == 2
        
        # Verify segmentation was created in database
        segmentation_id = response_data["segmentation_id"]
        segmentation = db_session.query(CellSegmentation).filter_by(id=segmentation_id).first()
        assert segmentation is not None
        assert segmentation.image_id == image.id
        assert segmentation.cell_count == 2

    def test_segment_invalid_image_integration(self, integration_client, timed_request):
        """Test segmentation with non-existent image returns 404."""
        fake_image_id = uuid4()
        
        response, benchmark = timed_request(
            "POST",
            f"/api/v1/imaging/{fake_image_id}/segment",
            "test_segment_invalid_image",
            json={"model": "cyto", "diameter": 30}
        )
        
        assert response.status_code == 404
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"

    @patch('amprenta_rag.api.routers.imaging.segment_batch.delay')
    def test_segment_batch_queues_task_integration(self, mock_task, integration_client, 
                                                 db_session, timed_request):
        """Test batch segmentation queues Celery task."""
        # Create real images in database
        images = []
        for i in range(3):
            image = MicroscopyImage(
                id=uuid4(),
                channel="DAPI",
                z_slice=0,
                timepoint=0,
                width=100,
                height=100,
                pixel_size_um=0.325,
                image_path=f"images/test/image_{i}.png",
                created_at=datetime.utcnow()
            )
            images.append(image)
        
        db_session.add_all(images)
        db_session.commit()
        
        # Mock Celery task (external service)
        mock_task.return_value.id = "task-123"
        
        image_ids = [str(img.id) for img in images]
        
        response, benchmark = timed_request(
            "POST",
            "/api/v1/imaging/segment/batch",
            "test_segment_batch",
            json={
                "image_ids": image_ids,
                "model": "cyto",
                "diameter": 30
            }
        )
        
        assert response.status_code == 202
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert "task_id" in response_data
        assert len(response_data["image_ids"]) == 3
        
        # Verify Celery task was called
        mock_task.assert_called_once()

    def test_get_image_metadata_integration(self, integration_client, db_session, timed_request):
        """Test retrieving image metadata from database."""
        # Create real image in database
        image = MicroscopyImage(
            id=uuid4(),
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=100,
            height=100,
            pixel_size_um=0.325,
            image_path="images/test/DAPI_z000_t000.png",
            created_at=datetime.utcnow()
        )
        db_session.add(image)
        db_session.commit()
        db_session.refresh(image)
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/imaging/{image.id}/metadata",
            "test_get_image_metadata"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert response_data["id"] == str(image.id)
        assert response_data["channel"] == "DAPI"
        assert response_data["width"] == 100
        assert response_data["height"] == 100
        assert response_data["pixel_size_um"] == 0.325

    def test_get_segmentation_integration(self, integration_client, db_session, timed_request):
        """Test retrieving segmentation data from database."""
        # Create real image and segmentation in database
        image = MicroscopyImage(
            id=uuid4(),
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=100,
            height=100,
            image_path="images/test/DAPI_z000_t000.png",
            created_at=datetime.utcnow()
        )
        db_session.add(image)
        db_session.commit()
        
        segmentation = CellSegmentation(
            id=uuid4(),
            image_id=image.id,
            model="cyto",
            diameter=30,
            cell_count=5,
            mask_path="masks/test_mask.npy",
            created_at=datetime.utcnow()
        )
        db_session.add(segmentation)
        db_session.commit()
        db_session.refresh(segmentation)
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/imaging/{image.id}/segmentation/{segmentation.id}",
            "test_get_segmentation"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert response_data["id"] == str(segmentation.id)
        assert response_data["cell_count"] == 5
        assert response_data["model"] == "cyto"

    @patch('amprenta_rag.api.routers.imaging.get_feature_extractor')
    def test_get_features_integration(self, mock_extractor, integration_client, 
                                    db_session, timed_request):
        """Test extracting features from real segmentation."""
        # Create real image, segmentation, and features in database
        image = MicroscopyImage(
            id=uuid4(),
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=100,
            height=100,
            image_path="images/test/DAPI_z000_t000.png",
            created_at=datetime.utcnow()
        )
        db_session.add(image)
        db_session.commit()
        
        segmentation = CellSegmentation(
            id=uuid4(),
            image_id=image.id,
            model="cyto",
            cell_count=2,
            mask_path="masks/test_mask.npy",
            created_at=datetime.utcnow()
        )
        db_session.add(segmentation)
        db_session.commit()
        
        # Create real features in database
        features = [
            CellFeature(
                id=uuid4(),
                segmentation_id=segmentation.id,
                cell_id=1,
                area=100.0,
                perimeter=40.0,
                circularity=0.8,
                centroid_x=50.0,
                centroid_y=60.0
            ),
            CellFeature(
                id=uuid4(),
                segmentation_id=segmentation.id,
                cell_id=2,
                area=120.0,
                perimeter=45.0,
                circularity=0.7,
                centroid_x=70.0,
                centroid_y=80.0
            )
        ]
        db_session.add_all(features)
        db_session.commit()
        
        # Mock feature extractor (external processing)
        mock_extractor.return_value.extract_morphology_features.return_value = []
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/imaging/{image.id}/segmentation/{segmentation.id}/features",
            "test_get_features"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert len(response_data["features"]) == 2
        assert response_data["features"][0]["area"] == 100.0
        assert response_data["features"][1]["area"] == 120.0

    @patch('amprenta_rag.api.routers.imaging.get_feature_extractor')
    def test_get_well_summary_integration(self, mock_extractor, integration_client,
                                        db_session, test_well, timed_request):
        """Test getting well summary with real database data."""
        # Create real images for the well
        images = []
        for channel in ["DAPI", "GFP"]:
            image = MicroscopyImage(
                id=uuid4(),
                well_id=test_well.id,
                channel=channel,
                z_slice=0,
                timepoint=0,
                width=100,
                height=100,
                image_path=f"images/test/{channel}_z000_t000.png",
                created_at=datetime.utcnow()
            )
            images.append(image)
        
        db_session.add_all(images)
        db_session.commit()
        
        # Mock feature extractor (external processing)
        mock_aggregated = MagicMock()
        mock_aggregated.to_dict.return_value = {
            "cell_count": 50,
            "morphology_stats": {"area": {"mean": 100.0, "std": 10.0}},
            "intensity_stats": {"DAPI": {"mean": 150.0}}
        }
        mock_extractor.return_value.aggregate_to_well.return_value = mock_aggregated
        mock_extractor.return_value.calculate_well_metrics.return_value = {
            "cell_count": 50, 
            "mean_cell_area": 100.0
        }
        
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/imaging/well/{test_well.id}/summary",
            "test_get_well_summary"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert response_data["well_id"] == str(test_well.id)
        assert response_data["image_count"] == 2
        assert len(response_data["channels"]) == 2
        assert "DAPI" in response_data["channels"]
        assert "GFP" in response_data["channels"]

    def test_get_well_summary_no_images_integration(self, integration_client, db_session, 
                                                   test_well, timed_request):
        """Test well summary with no images returns appropriate response."""
        response, benchmark = timed_request(
            "GET",
            f"/api/v1/imaging/well/{test_well.id}/summary",
            "test_get_well_summary_no_images"
        )
        
        assert response.status_code == 200
        assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms"
        
        response_data = response.json()
        assert response_data["well_id"] == str(test_well.id)
        assert response_data["image_count"] == 0
        assert response_data["channels"] == []
        assert response_data["features"] is None
