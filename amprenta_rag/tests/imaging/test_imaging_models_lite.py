"""Unit tests for imaging models and storage (SQLite compatible)."""

from __future__ import annotations

import json
import uuid
import numpy as np
import tempfile
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy import create_engine, Column, String, DateTime, UUID as SQLAlchemyUUID
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.postgresql import UUID

from amprenta_rag.database.base import Base
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.storage import ImageStorage, LocalStorageBackend, S3StorageBackend


def generate_uuid():
    return uuid.uuid4()


# Simplified models for testing (SQLite compatible)
class TestHTSCampaign(Base):
    """Simplified HTS campaign for testing."""
    __tablename__ = "test_hts_campaigns"
    
    id = Column(SQLAlchemyUUID(as_uuid=True), primary_key=True, default=generate_uuid)
    campaign_id = Column(String(200), nullable=False)
    campaign_name = Column(String(500), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


class TestHTSPlate(Base):
    """Simplified HTS plate for testing."""
    __tablename__ = "test_hts_plates"
    
    id = Column(SQLAlchemyUUID(as_uuid=True), primary_key=True, default=generate_uuid)
    campaign_id = Column(SQLAlchemyUUID(as_uuid=True), nullable=False)
    barcode = Column(String(200), nullable=False)
    plate_format = Column(String(50), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


class TestHTSWell(Base):
    """Simplified HTS well for testing."""
    __tablename__ = "test_hts_wells"
    
    id = Column(SQLAlchemyUUID(as_uuid=True), primary_key=True, default=generate_uuid)
    plate_id = Column(SQLAlchemyUUID(as_uuid=True), nullable=False)
    position = Column(String(10), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


@pytest.fixture
def db_session():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    
    # Create only the tables we need for testing
    TestHTSCampaign.__table__.create(engine)
    TestHTSPlate.__table__.create(engine)
    TestHTSWell.__table__.create(engine)
    MicroscopyImage.__table__.create(engine)
    CellSegmentation.__table__.create(engine)
    CellFeature.__table__.create(engine)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    
    yield session
    
    session.close()


@pytest.fixture
def sample_campaign(db_session):
    """Create a sample HTS campaign for testing."""
    campaign_id = uuid.uuid4()
    campaign = TestHTSCampaign(
        id=campaign_id,
        campaign_id="CAMP-001",
        campaign_name="Test Campaign"
    )
    db_session.add(campaign)
    db_session.commit()
    return campaign_id


@pytest.fixture
def sample_plate(db_session, sample_campaign):
    """Create a sample HTS plate for testing."""
    plate_id = uuid.uuid4()
    plate = TestHTSPlate(
        id=plate_id,
        campaign_id=sample_campaign,  # Now this is just the UUID
        barcode="PLATE-001",
        plate_format="384"
    )
    db_session.add(plate)
    db_session.commit()
    return plate_id


@pytest.fixture
def sample_well(db_session, sample_plate):
    """Create a sample HTS well for testing."""
    well_id = uuid.uuid4()
    well = TestHTSWell(
        id=well_id,
        plate_id=sample_plate,  # Now this is just the UUID
        position="A03"
    )
    db_session.add(well)
    db_session.commit()
    return well_id


class TestImagingModels:
    """Test imaging and microscopy models."""

    def test_microscopy_image_creation(self, db_session, sample_well):
        """Test creating microscopy image with well FK."""
        image = MicroscopyImage(
            well_id=sample_well,
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=2048,
            height=2048,
            bit_depth=16,
            pixel_size_um=0.325,
            image_path="images/well_123/DAPI_z000_t000.tiff",
            image_metadata={"exposure_ms": 100, "gain": 2.0}
        )
        
        db_session.add(image)
        db_session.commit()
        
        # Verify image was created
        assert image.id is not None
        assert image.channel == "DAPI"
        assert image.width == 2048
        assert image.height == 2048
        assert image.well_id == sample_well
        assert image.image_metadata["exposure_ms"] == 100

    def test_cell_segmentation_creation(self, db_session, sample_well):
        """Test creating cell segmentation with image FK."""
        # First create an image
        image = MicroscopyImage(
            well_id=sample_well,
            channel="GFP",
            width=1024,
            height=1024,
            bit_depth=16,
            image_path="images/test.tiff"
        )
        db_session.add(image)
        db_session.commit()
        
        # Create segmentation
        segmentation = CellSegmentation(
            image_id=image.id,
            model_name="cellpose",
            model_version="2.0",
            cell_count=157,
            mask_path="masks/seg_123_mask.npy",
            parameters={"diameter": 30, "flow_threshold": 0.4},
            confidence_score=0.92
        )
        
        db_session.add(segmentation)
        db_session.commit()
        
        # Verify segmentation was created
        assert segmentation.id is not None
        assert segmentation.model_name == "cellpose"
        assert segmentation.cell_count == 157
        assert segmentation.image_id == image.id
        assert segmentation.parameters["diameter"] == 30

    def test_cell_feature_creation(self, db_session, sample_well):
        """Test creating cell feature with segmentation FK."""
        # Create image and segmentation first
        image = MicroscopyImage(
            well_id=sample_well,
            channel="Brightfield",
            width=512,
            height=512,
            bit_depth=8,
            image_path="images/test.tiff"
        )
        db_session.add(image)
        db_session.commit()
        
        segmentation = CellSegmentation(
            image_id=image.id,
            model_name="custom_unet",
            cell_count=45,
            mask_path="masks/test_mask.npy"
        )
        db_session.add(segmentation)
        db_session.commit()
        
        # Create cell feature
        feature = CellFeature(
            segmentation_id=segmentation.id,
            cell_id=1,
            area=1250.5,
            perimeter=145.2,
            circularity=0.78,
            eccentricity=0.45,
            centroid_x=256.3,
            centroid_y=128.7,
            intensity_features={
                "DAPI": {"mean": 125.5, "std": 23.1},
                "GFP": {"mean": 89.2, "std": 15.4}
            },
            texture_features={"contrast": 0.65, "homogeneity": 0.82}
        )
        
        db_session.add(feature)
        db_session.commit()
        
        # Verify feature was created
        assert feature.id is not None
        assert feature.cell_id == 1
        assert feature.area == 1250.5
        assert feature.segmentation_id == segmentation.id
        assert feature.intensity_features["DAPI"]["mean"] == 125.5

    def test_well_channel_index(self, db_session, sample_well):
        """Test querying by well_id + channel using index."""
        # Create multiple images for the same well
        channels = ["DAPI", "GFP", "RFP"]
        images = []
        
        for channel in channels:
            image = MicroscopyImage(
                well_id=sample_well,
                channel=channel,
                width=1024,
                height=1024,
                bit_depth=16,
                image_path=f"images/{channel}.tiff"
            )
            images.append(image)
            db_session.add(image)
        
        db_session.commit()
        
        # Query by well_id and channel (should use index)
        dapi_images = db_session.query(MicroscopyImage).filter_by(
            well_id=sample_well,
            channel="DAPI"
        ).all()
        
        assert len(dapi_images) == 1
        assert dapi_images[0].channel == "DAPI"
        
        # Query all images for well
        well_images = db_session.query(MicroscopyImage).filter_by(
            well_id=sample_well
        ).all()
        
        assert len(well_images) == 3
        image_channels = {img.channel for img in well_images}
        assert image_channels == {"DAPI", "GFP", "RFP"}

    def test_imaging_relationships(self, db_session, sample_well):
        """Test relationships between imaging models."""
        # Create image
        image = MicroscopyImage(
            well_id=sample_well,
            channel="DAPI",
            width=1024,
            height=1024,
            bit_depth=16,
            image_path="images/test.tiff"
        )
        db_session.add(image)
        db_session.commit()
        
        # Create segmentation
        segmentation = CellSegmentation(
            image_id=image.id,
            model_name="cellpose",
            cell_count=10,
            mask_path="masks/test_mask.npy"
        )
        db_session.add(segmentation)
        db_session.commit()
        
        # Create multiple features
        for cell_id in range(1, 4):  # 3 cells
            feature = CellFeature(
                segmentation_id=segmentation.id,
                cell_id=cell_id,
                area=100.0 + cell_id * 10,
                perimeter=50.0 + cell_id * 5
            )
            db_session.add(feature)
        
        db_session.commit()
        
        # Test relationships
        # Image -> Segmentations
        assert len(image.segmentations) == 1
        assert image.segmentations[0].id == segmentation.id
        
        # Segmentation -> Features
        assert len(segmentation.features) == 3
        cell_ids = {f.cell_id for f in segmentation.features}
        assert cell_ids == {1, 2, 3}


class TestImageStorage:
    """Test image storage functionality."""

    def test_image_storage_save_load(self):
        """Test saving and loading image files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create storage with local backend
            storage = ImageStorage.create_local(temp_dir)
            
            # Create test image data
            image_data = np.random.randint(0, 255, (100, 100), dtype=np.uint8)
            
            # Save image
            path = storage.save_image(
                image_data=image_data,
                well_id="well_123",
                channel="DAPI",
                z_slice=0,
                timepoint=0,
                format="png"
            )
            
            # Verify path format
            expected_path = "images/well_123/DAPI_z000_t000.png"
            assert path == expected_path
            
            # Load image back
            loaded_image = storage.load_image(path)
            
            # Verify image data (allowing for compression artifacts)
            assert loaded_image.shape == image_data.shape
            assert loaded_image.dtype == image_data.dtype

    def test_mask_storage_npy(self):
        """Test saving and loading NPY mask files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            storage = ImageStorage.create_local(temp_dir)
            
            # Create test mask data
            mask_data = np.random.randint(0, 100, (200, 200), dtype=np.int32)
            
            # Save mask
            path = storage.save_mask(
                mask_data=mask_data,
                segmentation_id="seg_456"
            )
            
            # Verify path format
            expected_path = "masks/seg_456_mask.npy"
            assert path == expected_path
            
            # Load mask back
            loaded_mask = storage.load_mask(path)
            
            # Verify exact match (no compression for NPY)
            assert np.array_equal(loaded_mask, mask_data)
            assert loaded_mask.shape == mask_data.shape
            assert loaded_mask.dtype == mask_data.dtype

    def test_storage_backend_methods(self):
        """Test storage backend utility methods."""
        with tempfile.TemporaryDirectory() as temp_dir:
            storage = ImageStorage.create_local(temp_dir)
            
            # Test file that doesn't exist
            assert not storage.exists("nonexistent/path.png")
            
            # Save a file
            test_data = np.ones((50, 50), dtype=np.uint8) * 128
            path = storage.save_image(test_data, "test_well", "TEST", format="png")
            
            # Test file exists
            assert storage.exists(path)
            
            # Test metadata extraction
            metadata = storage.get_image_metadata(path)
            assert metadata["width"] == 50
            assert metadata["height"] == 50
            
            # Test file deletion
            assert storage.delete_image(path)
            assert not storage.exists(path)

    @patch('amprenta_rag.imaging.storage.HAS_S3', True)
    @patch('boto3.Session')
    def test_s3_storage_backend(self, mock_session):
        """Test S3 storage backend initialization."""
        # Mock S3 client
        mock_client = MagicMock()
        mock_session.return_value.client.return_value = mock_client
        
        # Create S3 storage backend
        backend = S3StorageBackend(
            bucket_name="test-bucket",
            prefix="imaging/",
            aws_access_key_id="test-key",
            aws_secret_access_key="test-secret"
        )
        
        # Verify initialization
        assert backend.bucket_name == "test-bucket"
        assert backend.prefix == "imaging/"
        
        # Test key generation
        s3_key = backend._get_s3_key("images/test.png")
        assert s3_key == "imaging/images/test.png"

    def test_local_storage_backend_direct(self):
        """Test LocalStorageBackend directly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            backend = LocalStorageBackend(temp_dir)
            
            # Test save and load
            test_data = b"test image data"
            path = backend.save_file(test_data, "test/image.png")
            
            # Verify file was saved
            assert backend.exists("test/image.png")
            
            # Load and verify
            loaded_data = backend.load_file("test/image.png")
            assert loaded_data == test_data
            
            # Test delete
            assert backend.delete_file("test/image.png")
            assert not backend.exists("test/image.png")
