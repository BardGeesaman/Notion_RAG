"""
Tests for biophysical assay API endpoints.

This test suite covers core biophysical API functionality using simplified mocking:
- SPR, MST, DSC file upload
- Basic endpoint validation

Note: Complex database interaction tests are deferred to integration tests 
due to SQLite threading limitations in FastAPI tests.
"""

from io import BytesIO
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.models.auth import User


@pytest.fixture
def mock_user():
    """Create a mock user for testing."""
    user = MagicMock(spec=User)
    user.id = uuid4()
    user.username = "testuser"
    user.email = "test@test.com"
    return user


@pytest.fixture
def client(mock_user):
    """Create test client with dependency overrides."""
    
    def override_get_db():
        return MagicMock()
    
    def override_get_current_user():
        return mock_user
    
    app.dependency_overrides[get_database_session] = override_get_db
    app.dependency_overrides[get_current_user] = override_get_current_user
    
    try:
        with TestClient(app) as client:
            yield client
    finally:
        app.dependency_overrides.clear()


# Test data fixtures
@pytest.fixture
def mock_spr_file():
    """Create a mock SPR CSV file for testing."""
    csv_content = "Time,Response,Cycle,Concentration\n0.0,0.0,1,1e-9\n1.0,10.0,1,1e-9\n2.0,20.0,1,1e-9\n3.0,15.0,1,1e-9\n4.0,5.0,1,1e-9\n"
    return BytesIO(csv_content.encode())


@pytest.fixture
def mock_mst_file():
    """Create a mock MST CSV file for testing."""
    csv_content = "Concentration,Fnorm,Fnorm_Error\n1e-12,1.0,0.1\n1e-11,0.95,0.1\n1e-10,0.8,0.1\n1e-9,0.5,0.1\n1e-8,0.2,0.1\n"
    return BytesIO(csv_content.encode())


@pytest.fixture
def mock_dsc_file():
    """Create a mock DSC CSV file for testing."""
    csv_content = "Temperature,Cp\n25.0,1.5\n30.0,1.6\n35.0,1.8\n40.0,2.2\n45.0,3.5\n50.0,5.2\n55.0,4.8\n60.0,3.1\n65.0,2.0\n70.0,1.7\n"
    return BytesIO(csv_content.encode())


# ============================================================================
# Core API Tests (3 tests covering main functionality)
# ============================================================================

def test_upload_spr_file_success(client, mock_spr_file):
    """Test successful SPR file upload."""
    
    with patch("amprenta_rag.api.routers.biophysical.ingest_spr_file") as mock_ingest:
        # Mock the ingest service to return a test SPR experiment
        mock_experiment = MagicMock()
        mock_experiment.id = uuid4()
        mock_experiment.experiment_name = "Test SPR"
        mock_experiment.processing_status = "pending"
        mock_ingest.return_value = mock_experiment
        
        response = client.post(
            "/api/v1/biophysical/spr/upload",
            files={"file": ("test_spr.csv", mock_spr_file, "text/csv")},
            data={"target_name": "Test Target"}
        )
    
    assert response.status_code == 201
    data = response.json()
    assert "experiment_id" in data
    assert data["filename"] == "test_spr.csv"
    assert data["processing_status"] == "pending"
    assert "Test SPR" in data["message"]


def test_upload_mst_file_success(client, mock_mst_file):
    """Test successful MST file upload."""
    
    with patch("amprenta_rag.api.routers.biophysical.ingest_mst_file") as mock_ingest:
        # Mock the ingest service to return a test MST experiment
        mock_experiment = MagicMock()
        mock_experiment.id = uuid4()
        mock_experiment.experiment_name = "Test MST"
        mock_experiment.processing_status = "pending"
        mock_ingest.return_value = mock_experiment
        
        response = client.post(
            "/api/v1/biophysical/mst/upload",
            files={"file": ("test_mst.csv", mock_mst_file, "text/csv")},
            data={"target_name": "Test Target"}
        )
    
    assert response.status_code == 201
    data = response.json()
    assert "experiment_id" in data
    assert data["filename"] == "test_mst.csv"
    assert data["processing_status"] == "pending"
    assert "Test MST" in data["message"]


def test_upload_dsc_file_success(client, mock_dsc_file):
    """Test successful DSC file upload."""
    
    with patch("amprenta_rag.api.routers.biophysical.ingest_dsc_file") as mock_ingest:
        # Mock the ingest service to return a test DSC experiment
        mock_experiment = MagicMock()
        mock_experiment.id = uuid4()
        mock_experiment.experiment_name = "Test DSC"
        mock_experiment.processing_status = "pending"
        mock_ingest.return_value = mock_experiment
        
        response = client.post(
            "/api/v1/biophysical/dsc/upload",
            files={"file": ("test_dsc.csv", mock_dsc_file, "text/csv")},
            data={"protein_name": "Test Protein"}
        )
    
    assert response.status_code == 201
    data = response.json()
    assert "experiment_id" in data
    assert data["filename"] == "test_dsc.csv"
    assert data["processing_status"] == "pending"
    assert "Test DSC" in data["message"]


# Note: Additional tests for listing, retrieval, refitting, and cross-assay comparison
# are deferred to integration tests in Batch 7 due to SQLite threading limitations
# in FastAPI TestClient. The upload tests above cover the core API functionality.
