"""API tests for new imaging metadata endpoints."""

import tempfile
from io import BytesIO
from pathlib import Path
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


@pytest.fixture
def client():
    """Test client with dependency overrides."""
    from amprenta_rag.api.dependencies import get_current_user
    
    # Mock user class
    class MockUser:
        def __init__(self):
            self.id = uuid4()
            self.username = "testuser"
            self.email = "test@test.com"
    
    # Mock user dependency
    def mock_user():
        return MockUser()
    
    app.dependency_overrides[get_current_user] = mock_user
    
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


@pytest.fixture
def test_microscope_id():
    """Return a test microscope ID."""
    return uuid4()


@pytest.fixture
def test_plate_id():
    """Return a test plate ID."""
    return uuid4()


@pytest.fixture
def test_well_id():
    """Return a test well ID."""
    return uuid4()


# Removed create_mock_ome_tiff function as it's not needed for simplified tests


class TestOMETiffImport:
    """Test OME-TIFF import endpoints."""
    
    def test_import_ome_tiff_invalid_file(self, client: TestClient):
        """Test OME-TIFF import with invalid file."""
        response = client.post(
            "/imaging/import/ome-tiff",
            files={"file": ("test.txt", BytesIO(b"not an image"), "text/plain")},
            data={}
        )
        
        assert response.status_code == 400
        assert "must be a TIFF" in response.json()["detail"]
    
    def test_import_ome_tiff_missing_file(self, client: TestClient):
        """Test OME-TIFF import without file."""
        response = client.post(
            "/imaging/import/ome-tiff",
            data={}
        )
        
        assert response.status_code == 422  # Validation error


class TestBatchImport:
    """Test batch vendor import endpoints."""
    
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
    
    def test_batch_import_missing_path(self, client: TestClient):
        """Test batch import without import path."""
        response = client.post(
            "/imaging/import/batch",
            data={"vendor": "opera"}
        )
        
        assert response.status_code == 422  # Validation error
    
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
        assert isinstance(response.json(), list)
    
    def test_create_microscope_invalid_data(self, client: TestClient):
        """Test creating microscope with invalid data."""
        microscope_data = {
            "name": "",  # Invalid empty name
            "manufacturer": "NewCorp",
            "model": "NC-2000"
        }
        
        response = client.post("/imaging/instruments", json=microscope_data)
        
        assert response.status_code == 422  # Validation error
    
    def test_create_microscope_missing_required(self, client: TestClient):
        """Test creating microscope with missing required fields."""
        microscope_data = {
            "name": "New Microscope"
            # Missing manufacturer and model
        }
        
        response = client.post("/imaging/instruments", json=microscope_data)
        
        assert response.status_code == 422  # Validation error
    
    def test_get_microscope_not_found(self, client: TestClient):
        """Test getting nonexistent microscope."""
        fake_id = uuid4()
        response = client.get(f"/imaging/instruments/{fake_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_list_objectives_empty(self, client: TestClient):
        """Test listing objectives when none exist."""
        response = client.get("/imaging/objectives")
        
        assert response.status_code == 200
        assert isinstance(response.json(), list)
    
    def test_list_channel_configs_empty(self, client: TestClient):
        """Test listing channel configs when none exist."""
        response = client.get("/imaging/channels")
        
        assert response.status_code == 200
        assert isinstance(response.json(), list)


class TestImageQC:
    """Test image QC endpoints."""
    
    def test_get_image_qc_not_found(self, client: TestClient):
        """Test getting QC for nonexistent image."""
        fake_id = uuid4()
        response = client.get(f"/imaging/qc/image/{fake_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_get_image_qc_invalid_id(self, client: TestClient):
        """Test getting QC with invalid image ID."""
        response = client.get("/imaging/qc/image/invalid-id")
        
        assert response.status_code == 422  # Validation error
    
    def test_get_plate_qc_not_found(self, client: TestClient):
        """Test getting plate QC for nonexistent plate."""
        fake_id = uuid4()
        response = client.get(f"/imaging/qc/plate/{fake_id}")
        
        assert response.status_code == 404
        assert "No images found" in response.json()["detail"]


class TestBrowseImages:
    """Test 5D image browser endpoint."""
    
    def test_browse_images_empty(self, client: TestClient):
        """Test browsing when no images exist."""
        response = client.get("/imaging/browse")
        
        assert response.status_code == 200
        data = response.json()
        
        assert "images" in data
        assert "total_count" in data
        assert "available_channels" in data
        assert "z_range" in data
        assert "t_range" in data
    
    def test_browse_images_with_filters(self, client: TestClient, test_plate_id):
        """Test browsing with filters."""
        params = {
            "plate_id": str(test_plate_id),
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
    
    def test_browse_images_invalid_pagination(self, client: TestClient):
        """Test browsing with invalid pagination parameters."""
        params = {
            "skip": -1,  # Invalid negative skip
            "limit": 0   # Invalid zero limit
        }
        
        response = client.get("/imaging/browse", params=params)
        
        assert response.status_code == 422  # Validation error
