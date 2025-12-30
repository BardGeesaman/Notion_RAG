"""Unit tests for imaging models and storage."""

from __future__ import annotations

import json
import uuid
import numpy as np
import tempfile
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from amprenta_rag.database.base import Base
from amprenta_rag.models.chemistry import HTSCampaign, HTSPlate, HTSWell, HTSResult, Compound
from amprenta_rag.imaging.models import MicroscopyImage, CellSegmentation, CellFeature
from amprenta_rag.imaging.storage import ImageStorage, LocalStorageBackend, S3StorageBackend


@pytest.fixture
def db_session():
    """Create in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:", echo=False)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    
    yield session
    
    session.close()


@pytest.fixture
def sample_compound(db_session):
    """Create a sample compound for testing."""
    compound = Compound(
        id=uuid.uuid4(),
        compound_id="COMP-001",
        smiles="CCO",  # Ethanol
        molecular_weight=46.07
    )
    db_session.add(compound)
    db_session.commit()
    db_session.expunge(compound)
    return compound


@pytest.fixture
def sample_campaign(db_session):
    """Create a sample HTS campaign for testing."""
    campaign = HTSCampaign(
        id=uuid.uuid4(),
        campaign_id="CAMP-001",
        campaign_name="Test Campaign",
        assay_type="Cell viability",
        target="CDK2"
    )
    db_session.add(campaign)
    db_session.commit()
    db_session.expunge(campaign)
    return campaign


@pytest.fixture
def sample_plate(db_session, sample_campaign):
    """Create a sample HTS plate for testing."""
    plate = HTSPlate(
        id=uuid.uuid4(),
        campaign_id=sample_campaign.id,
        barcode="PLATE-001",
        plate_format="384",
        layout_metadata={"control_wells": ["A01", "A02"]}
    )
    db_session.add(plate)
    db_session.commit()
    db_session.expunge(plate)
    return plate


@pytest.fixture
def sample_well(db_session, sample_plate, sample_compound):
    """Create a sample HTS well for testing."""
    well = HTSWell(
        id=uuid.uuid4(),
        plate_id=sample_plate.id,
        position="A03",
        compound_id=sample_compound.id,
        concentration=10.0,
        treatment_metadata={"treatment_time": 24}
    )
    db_session.add(well)
    db_session.commit()
    db_session.expunge(well)
    return well


class TestHTSPlateModels:
    """Test HTS plate hierarchy models."""

    def test_hts_plate_creation(self, db_session, sample_campaign):
        """Test creating HTS plate with campaign FK."""
        plate = HTSPlate(
            campaign_id=sample_campaign.id,
            barcode="TEST-PLATE-001",
            plate_format="1536",
            layout_metadata={"rows": 32, "cols": 48}
        )
        
        db_session.add(plate)
        db_session.commit()
        
        # Verify plate was created
        assert plate.id is not None
        assert plate.barcode == "TEST-PLATE-001"
        assert plate.plate_format == "1536"
        assert plate.layout_metadata["rows"] == 32
        assert plate.campaign_id == sample_campaign.id

    def test_hts_well_creation(self, db_session, sample_plate, sample_compound):
        """Test creating HTS well with plate FK and position."""
        well = HTSWell(
            plate_id=sample_plate.id,
            position="B12",
            compound_id=sample_compound.id,
            concentration=5.0,
            treatment_metadata={"dmso_percent": 0.1}
        )
        
        db_session.add(well)
        db_session.commit()
        
        # Verify well was created
        assert well.id is not None
        assert well.position == "B12"
        assert well.concentration == 5.0
        assert well.plate_id == sample_plate.id
        assert well.compound_id == sample_compound.id

    def test_hts_plate_well_relationship(self, db_session, sample_plate, sample_compound):
        """Test 1:N relationship between plate and wells."""
        # Create multiple wells for the plate
        wells = []
        for i, pos in enumerate(["A01", "A02", "A03"]):
            well = HTSWell(
                plate_id=sample_plate.id,
                position=pos,
                compound_id=sample_compound.id if i > 0 else None,  # A01 is control
                concentration=10.0 if i > 0 else None
            )
            wells.append(well)
            db_session.add(well)
        
        db_session.commit()
        
        # Query plate and verify wells relationship
        plate = db_session.query(HTSPlate).filter_by(id=sample_plate.id).first()
        assert len(plate.wells) == 3
        
        # Verify positions
        positions = {w.position for w in plate.wells}
        assert positions == {"A01", "A02", "A03"}

    def test_hts_result_well_fk(self, db_session, sample_campaign, sample_compound, sample_well):
        """Test optional well_id FK on HTSResult."""
        # Create HTS result with well FK
        result = HTSResult(
            result_id="RESULT-001",
            campaign_id=sample_campaign.id,
            compound_id=sample_compound.id,
            well_id=sample_well.id,  # P3 fix: well FK
            well_position=sample_well.position,
            raw_value=0.85,
            normalized_value=0.75,
            hit_flag=True
        )
        
        db_session.add(result)
        db_session.commit()
        
        # Verify result was created with well relationship
        assert result.id is not None
        assert result.well_id == sample_well.id
        assert result.well_position == "A03"
        
        # Test backward compatibility - result without well_id
        result_legacy = HTSResult(
            result_id="RESULT-002",
            campaign_id=sample_campaign.id,
            compound_id=sample_compound.id,
            well_id=None,  # Nullable for backward compat
            well_position="B01",
            raw_value=0.45
        )
        
        db_session.add(result_legacy)
        db_session.commit()
        
        assert result_legacy.id is not None
        assert result_legacy.well_id is None


class TestImagingModels:
    """Test imaging and microscopy models."""

    def test_microscopy_image_creation(self, db_session, sample_well):
        """Test creating microscopy image with well FK."""
        image = MicroscopyImage(
            well_id=sample_well.id,
            channel="DAPI",
            z_slice=0,
            timepoint=0,
            width=2048,
            height=2048,
            bit_depth=16,
            pixel_size_um=0.325,
            image_path="images/well_123/DAPI_z000_t000.tiff",
            metadata={"exposure_ms": 100, "gain": 2.0}
        )
        
        db_session.add(image)
        db_session.commit()
        
        # Verify image was created
        assert image.id is not None
        assert image.channel == "DAPI"
        assert image.width == 2048
        assert image.height == 2048
        assert image.well_id == sample_well.id
        assert image.metadata["exposure_ms"] == 100

    def test_cell_segmentation_creation(self, db_session, sample_well):
        """Test creating cell segmentation with image FK."""
        # First create an image
        image = MicroscopyImage(
            well_id=sample_well.id,
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
            well_id=sample_well.id,
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
                well_id=sample_well.id,
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
            well_id=sample_well.id,
            channel="DAPI"
        ).all()
        
        assert len(dapi_images) == 1
        assert dapi_images[0].channel == "DAPI"
        
        # Query all images for well
        well_images = db_session.query(MicroscopyImage).filter_by(
            well_id=sample_well.id
        ).all()
        
        assert len(well_images) == 3
        image_channels = {img.channel for img in well_images}
        assert image_channels == {"DAPI", "GFP", "RFP"}


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
