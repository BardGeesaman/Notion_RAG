"""API tests for new imaging metadata endpoints."""

import tempfile
from io import BytesIO
from pathlib import Path
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from amprenta_rag.api.main import app
from amprenta_rag.imaging.models_metadata import (
    Microscope, Objective, ChannelConfig, ImageFileSet
)
from amprenta_rag.models.chemistry import HTSWell, HTSPlate
from amprenta_rag.models.auth import User


@pytest.fixture
def client():
    """Test client with dependency overrides."""
    from amprenta_rag.api.dependencies import get_current_user
    
    # Mock user
    def mock_user():
        return User(id=uuid4(), username="testuser", email="test@test.com")
    
    app.dependency_overrides[get_current_user] = mock_user
    
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


@pytest.fixture
def test_microscope(db_session: Session):
    """Create test microscope."""
    microscope = Microscope(
        name="Test Microscope",
        manufacturer="Test Corp",
        model="TM-1000",
        serial_number="TM001",
        facility_location="Lab A"
    )
    db_session.add(microscope)
    db_session.commit()
    db_session.refresh(microscope)
    return microscope


@pytest.fixture
def test_plate(db_session: Session):
    """Create test plate."""
    plate = HTSPlate(
        plate_id="TEST_PLATE",
        plate_format=96,
        plate_type="imaging"
    )
    db_session.add(plate)
    db_session.commit()
    db_session.refresh(plate)
    return plate


@pytest.fixture
def test_well(db_session: Session, test_plate):
    """Create test well."""
    well = HTSWell(
        plate_id=test_plate.id,
        well_position="A01",
        row_index=0,
        col_index=0
    )
    db_session.add(well)
    db_session.commit()
    db_session.refresh(well)
    return well


def create_mock_ome_tiff():
    """Create a minimal mock OME-TIFF file."""
    import numpy as np
    import tifffile
    
    # Create small test image
    data = np.random.randint(0, 65535, (64, 64), dtype=np.uint16)
    
    # Create OME-XML metadata
    metadata = {
        'axes': 'YX',
        'PhysicalSizeX': 0.325,
        'PhysicalSizeY': 0.325,
    }
    
    # Write to bytes buffer
    buffer = BytesIO()
    tifffile.imwrite(buffer, data, ome=True, metadata=metadata)
    buffer.seek(0)
    
    return buffer.getvalue()


class TestOMETiffImport:
    """Test OME-TIFF import endpoints."""
    
    def test_import_ome_tiff_success(self, client: TestClient, test_well):
        """Test successful OME-TIFF import."""
        # Create mock OME-TIFF file
        tiff_data = create_mock_ome_tiff()
        
        response = client.post(
            "/imaging/import/ome-tiff",
            files={"file": ("test.ome.tiff", BytesIO(tiff_data), "image/tiff")},
            data={"well_id": str(test_well.id)}
        )
        
        assert response.status_code == 201
        data = response.json()
        
        assert "image_id" in data
        assert data["filename"] == "test.ome.tiff"
        assert "dimensions" in data
        assert "channels" in data
        assert "ome_metadata" in data
    
    def test_import_ome_tiff_invalid_file(self, client: TestClient):
        """Test OME-TIFF import with invalid file."""
        response = client.post(
            "/imaging/import/ome-tiff",
            files={"file": ("test.txt", BytesIO(b"not an image"), "text/plain")},
            data={}
        )
        
        assert response.status_code == 400
        assert "must be a TIFF" in response.json()["detail"]
    
    def test_import_ome_tiff_standalone(self, client: TestClient, test_plate):
        """Test standalone OME-TIFF import without well."""
        tiff_data = create_mock_ome_tiff()
        
        response = client.post(
            "/imaging/import/ome-tiff",
            files={"file": ("standalone.tiff", BytesIO(tiff_data), "image/tiff")},
            data={"plate_id": str(test_plate.id)}
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["filename"] == "standalone.tiff"


class TestBatchImport:
    """Test batch vendor import endpoints."""
    
    def test_batch_import_request(self, client: TestClient, test_plate):
        """Test batch import request."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create mock vendor export directory
            export_dir = Path(tmp_dir) / "opera_export"
            export_dir.mkdir()
            
            # Create minimal Index.xml
            index_xml = """<?xml version="1.0" encoding="UTF-8"?>
<OperaExport>
    <PlateBarcode>TEST_PLATE</PlateBarcode>
    <PlateType>96</PlateType>
</OperaExport>"""
            (export_dir / "Index.xml").write_text(index_xml)
            
            # Create Images directory
            images_dir = export_dir / "Images"
            images_dir.mkdir()
            
            response = client.post(
                "/imaging/import/batch",
                data={
                    "import_path": str(export_dir),
                    "plate_id": str(test_plate.id),
                    "vendor": "opera"
                }
            )
            
            assert response.status_code == 202
            data = response.json()
            
            assert "fileset_id" in data
            assert data["vendor"] == "opera"
            assert data["status"] in ["pending", "importing"]
    
    def test_batch_import_invalid_path(self, client: TestClient):
        """Test batch import with invalid path."""
        response = client.post(
            "/imaging/import/batch",
            data={
                "import_path": "/nonexistent/path",
                "vendor": "opera"
            }
        )
        
        assert response.status_code == 400
        assert "Failed to parse vendor export" in response.json()["detail"]
    
    def test_get_import_status(self, client: TestClient, db_session: Session, test_plate):
        """Test getting import job status."""
        # Create test fileset
        fileset = ImageFileSet(
            plate_id=test_plate.id,
            vendor="opera",
            import_path="/test/path",
            file_count=10,
            image_count=5,
            import_status="importing"
        )
        db_session.add(fileset)
        db_session.commit()
        db_session.refresh(fileset)
        
        response = client.get(f"/imaging/import/{fileset.id}/status")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["fileset_id"] == str(fileset.id)
        assert data["status"] == "importing"
        assert data["progress_percent"] == 50.0
        assert data["total_images"] == 10
        assert data["imported_count"] == 5
    
    def test_get_import_status_not_found(self, client: TestClient):
        """Test getting status for nonexistent import job."""
        fake_id = uuid4()
        response = client.get(f"/imaging/import/{fake_id}/status")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]


class TestInstrumentManagement:
    """Test instrument management endpoints."""
    
    def test_list_microscopes_empty(self, client: TestClient):
        """Test listing microscopes when none exist."""
        response = client.get("/imaging/instruments")
        
        assert response.status_code == 200
        assert response.json() == []
    
    def test_list_microscopes(self, client: TestClient, test_microscope):
        """Test listing microscopes."""
        response = client.get("/imaging/instruments")
        
        assert response.status_code == 200
        data = response.json()
        
        assert len(data) == 1
        microscope = data[0]
        assert microscope["name"] == "Test Microscope"
        assert microscope["manufacturer"] == "Test Corp"
        assert microscope["model"] == "TM-1000"
        assert microscope["is_active"] is True
    
    def test_create_microscope(self, client: TestClient):
        """Test creating a new microscope."""
        microscope_data = {
            "name": "New Microscope",
            "manufacturer": "NewCorp",
            "model": "NC-2000",
            "serial_number": "NC002",
            "facility_location": "Lab B"
        }
        
        response = client.post("/imaging/instruments", json=microscope_data)
        
        assert response.status_code == 201
        data = response.json()
        
        assert data["name"] == "New Microscope"
        assert data["manufacturer"] == "NewCorp"
        assert data["model"] == "NC-2000"
        assert data["serial_number"] == "NC002"
        assert data["is_active"] is True
    
    def test_create_microscope_duplicate_serial(self, client: TestClient, test_microscope):
        """Test creating microscope with duplicate serial number."""
        microscope_data = {
            "name": "Duplicate Microscope",
            "manufacturer": "DupeCorp",
            "model": "DUP-1000",
            "serial_number": test_microscope.serial_number  # Duplicate
        }
        
        response = client.post("/imaging/instruments", json=microscope_data)
        
        assert response.status_code == 400
        assert "already exists" in response.json()["detail"]
    
    def test_get_microscope(self, client: TestClient, test_microscope):
        """Test getting specific microscope details."""
        response = client.get(f"/imaging/instruments/{test_microscope.id}")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["id"] == str(test_microscope.id)
        assert data["name"] == test_microscope.name
        assert "objectives" in data
        assert "channels" in data
    
    def test_get_microscope_not_found(self, client: TestClient):
        """Test getting nonexistent microscope."""
        fake_id = uuid4()
        response = client.get(f"/imaging/instruments/{fake_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_list_objectives(self, client: TestClient, db_session: Session, test_microscope):
        """Test listing objectives."""
        # Create test objective
        objective = Objective(
            microscope_id=test_microscope.id,
            name="Test Objective 20x",
            magnification=20.0,
            numerical_aperture=0.75,
            immersion="air"
        )
        db_session.add(objective)
        db_session.commit()
        
        response = client.get("/imaging/objectives")
        
        assert response.status_code == 200
        data = response.json()
        
        assert len(data) == 1
        obj = data[0]
        assert obj["name"] == "Test Objective 20x"
        assert obj["magnification"] == 20.0
        assert obj["numerical_aperture"] == 0.75
    
    def test_list_channel_configs(self, client: TestClient, db_session: Session, test_microscope):
        """Test listing channel configurations."""
        # Create test channel config
        channel = ChannelConfig(
            microscope_id=test_microscope.id,
            channel_name="DAPI",
            fluorophore="Hoechst 33342",
            default_exposure_ms=100.0,
            default_gain=1.5
        )
        db_session.add(channel)
        db_session.commit()
        
        response = client.get("/imaging/channels")
        
        assert response.status_code == 200
        data = response.json()
        
        assert len(data) == 1
        ch = data[0]
        assert ch["channel_name"] == "DAPI"
        assert ch["fluorophore"] == "Hoechst 33342"
        assert ch["default_exposure_ms"] == 100.0


class TestImageQC:
    """Test image QC endpoints."""
    
    def test_get_image_qc_not_found(self, client: TestClient):
        """Test getting QC for nonexistent image."""
        fake_id = uuid4()
        response = client.get(f"/imaging/qc/image/{fake_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_get_plate_qc_no_images(self, client: TestClient, test_plate):
        """Test getting plate QC when no images exist."""
        response = client.get(f"/imaging/qc/plate/{test_plate.id}")
        
        assert response.status_code == 404
        assert "No images found" in response.json()["detail"]


class TestBrowseImages:
    """Test 5D image browser endpoint."""
    
    def test_browse_images_empty(self, client: TestClient):
        """Test browsing when no images exist."""
        response = client.get("/imaging/browse")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["images"] == []
        assert data["total_count"] == 0
        assert "available_channels" in data
        assert "z_range" in data
        assert "t_range" in data
    
    def test_browse_images_with_filters(self, client: TestClient, test_plate):
        """Test browsing with filters."""
        params = {
            "plate_id": str(test_plate.id),
            "channel": "DAPI",
            "limit": 50
        }
        
        response = client.get("/imaging/browse", params=params)
        
        assert response.status_code == 200
        data = response.json()
        
        assert "images" in data
        assert "total_count" in data
    
    def test_browse_images_pagination(self, client: TestClient):
        """Test browsing with pagination."""
        params = {
            "skip": 10,
            "limit": 20
        }
        
        response = client.get("/imaging/browse", params=params)
        
        assert response.status_code == 200
        data = response.json()
        
        assert "images" in data
        assert "total_count" in data
