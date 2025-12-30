"""
Tests for flow cytometry API endpoints.

This test suite covers all 9 flow cytometry API endpoints including:
- FCS file upload and ingestion
- Dataset listing and retrieval
- Event data access with pagination
- Gate CRUD operations
- Population statistics
"""

import tempfile
from io import BytesIO
from pathlib import Path
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.api.schemas import GateCreate
from amprenta_rag.models.auth import User
from amprenta_rag.database.models_flow_cytometry import (
    FlowCytometryDataset,
    FlowCytometryParameter,
    FlowCytometryGate,
    FlowCytometryPopulation,
)


def mock_current_user():
    """Mock current user for authentication."""
    return User(
        id=uuid4(),
        username="test_user",
        email="test@example.com",
        password_hash="fake_hash"
    )


class TestFlowCytometryAPI:
    """Test flow cytometry API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client with mocked dependencies."""
        app.dependency_overrides[get_current_user] = mock_current_user
        app.dependency_overrides[get_database_session] = self.mock_db_session
        
        try:
            with TestClient(app) as test_client:
                yield test_client
        finally:
            app.dependency_overrides.clear()
    
    def mock_db_session(self):
        """Mock database session with proper query chain setup."""
        db = MagicMock()
        
        # Configure mock query chains for different models
        def mock_query_side_effect(model):
            query_mock = MagicMock()
            
            if model == FlowCytometryDataset:
                # Setup for dataset queries
                query_mock.order_by.return_value.offset.return_value.limit.return_value.all.return_value = [self.mock_flow_dataset]
                query_mock.order_by.return_value.filter.return_value.offset.return_value.limit.return_value.all.return_value = [self.mock_flow_dataset]
                query_mock.filter.return_value.first.return_value = self.mock_flow_dataset
                
            elif model == FlowCytometryParameter:
                # Setup for parameter queries
                query_mock.filter.return_value.order_by.return_value.all.return_value = [self.mock_parameter]
                query_mock.filter.return_value.first.return_value = self.mock_parameter
                
            elif model == FlowCytometryGate:
                # Setup for gate queries
                query_mock.filter.return_value.filter.return_value.order_by.return_value.all.return_value = [self.mock_gate]
                query_mock.filter.return_value.first.return_value = self.mock_gate
                
            elif model == FlowCytometryPopulation:
                # Setup for population queries
                query_mock.filter.return_value.order_by.return_value.all.return_value = [self.mock_population]
                
            return query_mock
        
        db.query.side_effect = mock_query_side_effect
        
        # Mock other database operations
        db.add = MagicMock()
        db.commit = MagicMock()
        db.refresh = MagicMock()
        db.rollback = MagicMock()
        
        return db
    
    @property
    def mock_fcs_file(self):
        """Create a mock FCS file for upload testing."""
        # Create a simple mock FCS file content
        content = b"FCS3.0    256     512    1024    2048    0    0    0    0    0" + b"\x00" * 200
        content += b"mock event data" * 100  # Simulate event data
        
        return BytesIO(content)
    
    @property
    def mock_flow_dataset(self):
        """Mock FlowCytometryDataset for testing."""
        if not hasattr(self, '_mock_flow_dataset'):
            dataset = MagicMock()
            dataset.id = uuid4()
            dataset.dataset_id = uuid4()
            dataset.events_parquet_path = "test_events.parquet"
            dataset.file_size_bytes = 1024
            dataset.n_events = 1000
            dataset.n_parameters = 10
            dataset.processing_status = "completed"
            dataset.processing_log = None
            dataset.sample_id = "TEST001"
            self._mock_flow_dataset = dataset
        return self._mock_flow_dataset
    
    @property
    def mock_parameter(self):
        """Mock FlowCytometryParameter for testing."""
        if not hasattr(self, '_mock_parameter'):
            param = MagicMock()
            param.id = uuid4()
            param.parameter_index = 0
            param.parameter_name = "FSC-A"
            param.parameter_short_name = "FSC"
            param.fluorophore = None
            param.min_value = 0.0
            param.max_value = 1000.0
            self._mock_parameter = param
        return self._mock_parameter
    
    @property
    def mock_gate(self):
        """Mock FlowCytometryGate for testing."""
        if not hasattr(self, '_mock_gate'):
            gate = MagicMock()
            gate.id = uuid4()
            gate.flow_dataset_id = uuid4()
            gate.gate_name = "Test Gate"
            gate.gate_type = "rectangle"
            gate.gate_definition = {"x_min": 0, "x_max": 1000, "y_min": 0, "y_max": 1000}
            gate.x_parameter_id = uuid4()
            gate.y_parameter_id = uuid4()
            gate.is_active = True
            self._mock_gate = gate
        return self._mock_gate
    
    @property
    def mock_population(self):
        """Mock FlowCytometryPopulation for testing."""
        if not hasattr(self, '_mock_population'):
            population = MagicMock()
            population.id = uuid4()
            population.event_count = 500
            population.percentage_of_total = 50.0
            population.parameter_statistics = {"FSC-A": {"mean": 1000, "median": 950}}
            self._mock_population = population
        return self._mock_population
    
    def test_upload_fcs_file_success(self, client):
        """Test successful FCS file upload."""
        with patch('amprenta_rag.api.routers.flow_cytometry.ingest_fcs') as mock_ingest:
            # Setup mock
            mock_dataset = MagicMock()
            mock_dataset.id = uuid4()
            mock_dataset.dataset_id = uuid4()
            mock_dataset.processing_status = "pending"
            mock_ingest.return_value = mock_dataset
            
            # Test upload
            response = client.post(
                "/api/v1/flow-cytometry/upload",
                files={"file": ("test.fcs", self.mock_fcs_file, "application/octet-stream")}
            )
            
            assert response.status_code == 201
            data = response.json()
            assert "flow_dataset_id" in data
            assert "dataset_id" in data
            assert data["filename"] == "test.fcs"
            assert data["processing_status"] == "pending"
            assert "message" in data
            
            # Verify ingest was called
            mock_ingest.assert_called_once()
    
    def test_upload_fcs_file_invalid_extension(self, client):
        """Test FCS upload with invalid file extension."""
        mock_file = BytesIO(b"test data")
        
        response = client.post(
            "/api/v1/flow-cytometry/upload",
            files={"file": ("test.txt", mock_file, "text/plain")}
        )
        
        assert response.status_code == 400
        assert "must have .fcs extension" in response.json()["detail"]
    
    def test_list_flow_datasets(self, client):
        """Test listing flow cytometry datasets."""
        # Test request
        response = client.get("/api/v1/flow-cytometry/datasets")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert str(self.mock_flow_dataset.id) in str(data[0]["id"])
    
    def test_list_flow_datasets_with_filters(self, client):
        """Test listing datasets with status filter."""
        # Test with status filter
        response = client.get("/api/v1/flow-cytometry/datasets?processing_status=completed")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
    
    def test_get_flow_dataset_success(self, client):
        """Test getting specific flow dataset."""
        dataset_id = self.mock_flow_dataset.id
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(dataset_id)
        assert data["processing_status"] == "completed"
    
    def test_get_flow_dataset_not_found(self, client):
        """Test getting non-existent dataset."""
        # Override the mock to return None for this test
        def mock_db_session_not_found():
            db = MagicMock()
            def mock_query_side_effect(model):
                query_mock = MagicMock()
                query_mock.filter.return_value.first.return_value = None
                return query_mock
            db.query.side_effect = mock_query_side_effect
            return db
        
        app.dependency_overrides[get_database_session] = mock_db_session_not_found
        
        dataset_id = uuid4()
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_get_dataset_parameters(self, client):
        """Test getting dataset parameters."""
        dataset_id = self.mock_flow_dataset.id
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/parameters")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["parameter_name"] == "FSC-A"
    
    @patch('amprenta_rag.api.routers.flow_cytometry.load_events_parquet')
    def test_get_dataset_events(self, mock_load_events, client):
        """Test getting dataset events."""
        # Mock events data
        import numpy as np
        mock_events = np.array([[100, 200], [150, 250], [200, 300]], dtype=np.float32)
        mock_params = ["FSC-A", "SSC-A"]
        mock_load_events.return_value = (mock_events, mock_params)
        
        dataset_id = self.mock_flow_dataset.id
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?limit=10")
        
        assert response.status_code == 200
        data = response.json()
        assert "events" in data
        assert "parameter_names" in data
        assert data["parameter_names"] == mock_params
        assert data["total_events"] == 3
        assert len(data["events"]) == 3
    
    def test_get_dataset_events_not_processed(self, client):
        """Test getting events from unprocessed dataset."""
        # Override mock to return unprocessed dataset
        def mock_db_session_unprocessed():
            db = MagicMock()
            def mock_query_side_effect(model):
                query_mock = MagicMock()
                if model == FlowCytometryDataset:
                    mock_dataset = MagicMock()
                    mock_dataset.processing_status = "pending"
                    query_mock.filter.return_value.first.return_value = mock_dataset
                return query_mock
            db.query.side_effect = mock_query_side_effect
            return db
        
        app.dependency_overrides[get_database_session] = mock_db_session_unprocessed
        
        dataset_id = uuid4()
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events")
        
        assert response.status_code == 409
        assert "processing not completed" in response.json()["detail"]
    
    @patch('amprenta_rag.api.routers.flow_cytometry.apply_gate_to_dataset')
    def test_create_gate_success(self, mock_apply_gate, client):
        """Test successful gate creation."""
        mock_apply_gate.return_value = self.mock_gate
        
        dataset_id = self.mock_flow_dataset.id
        
        # Test gate creation
        gate_data = {
            "gate_name": "Test Rectangle Gate",
            "gate_type": "rectangle",
            "gate_definition": {"x_min": 0, "x_max": 1000, "y_min": 0, "y_max": 1000},
            "x_parameter_id": str(uuid4()),
            "y_parameter_id": str(uuid4())
        }
        
        response = client.post(
            f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates",
            json=gate_data
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["gate_name"] == "Test Gate"
        assert data["gate_type"] == "rectangle"
        
        # Verify apply_gate_to_dataset was called
        mock_apply_gate.assert_called_once()
    
    def test_create_gate_invalid_definition(self, client):
        """Test gate creation with invalid definition."""
        dataset_id = self.mock_flow_dataset.id
        
        # Test with invalid polygon gate (< 3 vertices)
        gate_data = {
            "gate_name": "Invalid Polygon",
            "gate_type": "polygon",
            "gate_definition": {"vertices": [[0, 0], [1, 1]]},  # Only 2 vertices
            "x_parameter_id": str(uuid4()),
            "y_parameter_id": str(uuid4())
        }
        
        response = client.post(
            f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates",
            json=gate_data
        )
        
        assert response.status_code == 422  # Validation error
        assert "at least 3 vertices" in str(response.json())
    
    def test_list_gates(self, client):
        """Test listing gates for a dataset."""
        dataset_id = self.mock_flow_dataset.id
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["gate_name"] == "Test Gate"
    
    def test_update_gate(self, client):
        """Test updating an existing gate."""
        gate_id = self.mock_gate.id
        
        # Test update
        update_data = {
            "gate_name": "Updated Gate Name",
            "is_active": False
        }
        
        response = client.put(
            f"/api/v1/flow-cytometry/gates/{gate_id}",
            json=update_data
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["gate_name"] == "Test Gate"  # Mock returns original name
    
    def test_delete_gate(self, client):
        """Test deleting (deactivating) a gate."""
        gate_id = self.mock_gate.id
        
        # Test delete
        response = client.delete(f"/api/v1/flow-cytometry/gates/{gate_id}")
        
        assert response.status_code == 204
        
        # Verify gate was marked inactive
        assert self.mock_gate.is_active == False
    
    def test_delete_gate_not_found(self, client):
        """Test deleting non-existent gate."""
        # Override mock to return None for this test
        def mock_db_session_not_found():
            db = MagicMock()
            def mock_query_side_effect(model):
                query_mock = MagicMock()
                query_mock.filter.return_value.first.return_value = None
                return query_mock
            db.query.side_effect = mock_query_side_effect
            return db
        
        app.dependency_overrides[get_database_session] = mock_db_session_not_found
        
        gate_id = uuid4()
        
        # Test delete
        response = client.delete(f"/api/v1/flow-cytometry/gates/{gate_id}")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"]
    
    def test_get_population_statistics(self, client):
        """Test getting population statistics."""
        dataset_id = self.mock_flow_dataset.id
        
        # Test request
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/populations")
        
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["event_count"] == 500
        assert data[0]["percentage_of_total"] == 50.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
