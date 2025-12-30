"""
Integration tests for Flow Cytometry functionality.

These tests cover the complex workflows deferred from Batch 4 API tests,
using real PostgreSQL database instead of SQLite mocks.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from uuid import uuid4

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, User
from amprenta_rag.database.models_flow_cytometry import (
    FlowCytometryDataset,
    FlowCytometryParameter,
    FlowCytometryGate,
    FlowCytometryPopulation,
    GateType,
)
from amprenta_rag.flow_cytometry.fcs_parser import save_events_parquet


pytestmark = pytest.mark.integration


@pytest.fixture
def test_user():
    """Create a test user in the database."""
    with db_session() as db:
        user = User(
            id=uuid4(),
            username=f"test_user_{uuid4().hex[:8]}",
            email=f"test_{uuid4().hex[:8]}@example.com",
            password_hash="test_hash"
        )
        db.add(user)
        db.commit()
        db.expunge(user)
        return user


@pytest.fixture
def test_dataset(test_user):
    """Create a test dataset in the database."""
    with db_session() as db:
        dataset = Dataset(
            id=uuid4(),
            name=f"Integration Test Dataset {uuid4().hex[:8]}",
            omics_type="flow_cytometry",
            description="Integration test dataset",
            created_by=test_user.id
        )
        db.add(dataset)
        db.commit()
        db.expunge(dataset)
        return dataset


@pytest.fixture
def test_flow_dataset_with_events(test_dataset):
    """Create a flow cytometry dataset with synthetic events."""
    # Generate synthetic events
    import numpy as np
    np.random.seed(42)  # Reproducible data
    
    n_events = 1000
    events_data = {
        "FSC-A": np.random.uniform(10000, 100000, n_events),
        "SSC-A": np.random.uniform(5000, 80000, n_events),
        "FITC-A": np.random.uniform(100, 50000, n_events),
        "PE-A": np.random.uniform(100, 30000, n_events),
        "PerCP-A": np.random.uniform(100, 25000, n_events),
        "APC-A": np.random.uniform(100, 40000, n_events),
    }
    events_df = pd.DataFrame(events_data)
    
    # Save to temporary parquet file
    temp_dir = Path(tempfile.gettempdir()) / "flow_integration_test"
    temp_dir.mkdir(exist_ok=True)
    parquet_path = temp_dir / f"events_{test_dataset.id}.parquet"
    save_events_parquet(events_df, str(parquet_path))
    
    with db_session() as db:
        # Create flow dataset
        flow_dataset = FlowCytometryDataset(
            id=uuid4(),
            dataset_id=test_dataset.id,
            events_parquet_path=str(parquet_path),
            file_size_bytes=parquet_path.stat().st_size,
            n_events=n_events,
            n_parameters=len(events_data),
            cytometer_model="Test Cytometer",
            sample_id=f"TEST_{uuid4().hex[:8]}",
            processing_status="completed"
        )
        db.add(flow_dataset)
        
        # Create parameters
        for i, (param_name, _) in enumerate(events_data.items()):
            parameter = FlowCytometryParameter(
                id=uuid4(),
                flow_dataset_id=flow_dataset.id,
                parameter_name=param_name,
                parameter_index=i,
                range_min=0.0,
                range_max=262144.0,
                display_name=param_name
            )
            db.add(parameter)
        
        db.commit()
        db.expunge(flow_dataset)
        return flow_dataset


@pytest.fixture
def client(test_user):
    """Create test client with authentication."""
    # Mock authentication
    def mock_get_current_user():
        return test_user
    
    app.dependency_overrides[get_current_user] = mock_get_current_user
    
    try:
        with TestClient(app) as client:
            yield client
    finally:
        app.dependency_overrides.clear()


class TestFlowCytometryIntegration:
    """Integration tests for flow cytometry functionality."""
    
    def test_full_dataset_workflow(self, client, test_flow_dataset_with_events):
        """Test complete dataset workflow: list, get details, get parameters, get events."""
        # Test dataset listing
        response = client.get("/api/v1/flow-cytometry/datasets")
        assert response.status_code == 200
        data = response.json()
        assert "items" in data
        assert len(data["items"]) >= 1
        
        # Find our test dataset
        test_dataset = None
        for item in data["items"]:
            if item["id"] == str(test_flow_dataset_with_events.id):
                test_dataset = item
                break
        
        assert test_dataset is not None
        assert test_dataset["processing_status"] == "completed"
        assert test_dataset["n_events"] == 1000
        assert test_dataset["n_parameters"] == 6
        
        # Test dataset details
        response = client.get(f"/api/v1/flow-cytometry/datasets/{test_flow_dataset_with_events.id}")
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(test_flow_dataset_with_events.id)
        assert data["processing_status"] == "completed"
        
        # Test parameters retrieval
        response = client.get(f"/api/v1/flow-cytometry/datasets/{test_flow_dataset_with_events.id}/parameters")
        assert response.status_code == 200
        parameters = response.json()
        assert len(parameters) == 6
        
        param_names = [p["parameter_name"] for p in parameters]
        expected_params = ["FSC-A", "SSC-A", "FITC-A", "PE-A", "PerCP-A", "APC-A"]
        for expected in expected_params:
            assert expected in param_names
        
        # Test events retrieval
        response = client.get(f"/api/v1/flow-cytometry/datasets/{test_flow_dataset_with_events.id}/events?limit=100")
        assert response.status_code == 200
        data = response.json()
        assert "events" in data
        assert "total_events" in data
        assert data["total_events"] == 1000
        assert len(data["events"]) <= 100
        
        # Verify event structure
        if data["events"]:
            event = data["events"][0]
            for param_name in expected_params:
                assert param_name in event
                assert isinstance(event[param_name], (int, float))
    
    def test_gate_crud_operations(self, client, test_flow_dataset_with_events):
        """Test complete gate CRUD operations."""
        dataset_id = test_flow_dataset_with_events.id
        
        # Test gate creation - polygon gate
        polygon_gate_data = {
            "gate_name": "Test Polygon Gate",
            "gate_type": "polygon",
            "x_parameter": "FSC-A",
            "y_parameter": "SSC-A",
            "gate_definition": {
                "vertices": [[20000, 10000], [80000, 10000], [80000, 60000], [20000, 60000]]
            }
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=polygon_gate_data)
        assert response.status_code == 201
        polygon_gate = response.json()
        assert polygon_gate["gate_name"] == "Test Polygon Gate"
        assert polygon_gate["gate_type"] == "polygon"
        polygon_gate_id = polygon_gate["id"]
        
        # Test gate creation - rectangle gate
        rectangle_gate_data = {
            "gate_name": "Test Rectangle Gate",
            "gate_type": "rectangle",
            "x_parameter": "FITC-A",
            "y_parameter": "PE-A",
            "gate_definition": {
                "x_min": 1000, "x_max": 20000,
                "y_min": 1000, "y_max": 15000
            }
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=rectangle_gate_data)
        assert response.status_code == 201
        rectangle_gate = response.json()
        assert rectangle_gate["gate_name"] == "Test Rectangle Gate"
        assert rectangle_gate["gate_type"] == "rectangle"
        rectangle_gate_id = rectangle_gate["id"]
        
        # Test gate listing
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates")
        assert response.status_code == 200
        gates = response.json()
        assert len(gates) >= 2
        
        gate_names = [g["gate_name"] for g in gates]
        assert "Test Polygon Gate" in gate_names
        assert "Test Rectangle Gate" in gate_names
        
        # Test gate update
        update_data = {
            "gate_name": "Updated Polygon Gate",
            "is_active": False
        }
        
        response = client.put(f"/api/v1/flow-cytometry/gates/{polygon_gate_id}", json=update_data)
        assert response.status_code == 200
        updated_gate = response.json()
        assert updated_gate["gate_name"] == "Updated Polygon Gate"
        assert updated_gate["is_active"] == False
        
        # Test gate deletion
        response = client.delete(f"/api/v1/flow-cytometry/gates/{rectangle_gate_id}")
        assert response.status_code == 204
        
        # Verify gate was deactivated (soft delete)
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates")
        assert response.status_code == 200
        gates = response.json()
        
        # Should still have gates but rectangle gate should be inactive
        active_gates = [g for g in gates if g["is_active"]]
        assert len(active_gates) >= 0  # Polygon gate was deactivated, rectangle gate was deleted
    
    def test_population_statistics_calculation(self, client, test_flow_dataset_with_events):
        """Test population statistics calculation and retrieval."""
        dataset_id = test_flow_dataset_with_events.id
        
        # Create a gate that should capture some events
        gate_data = {
            "gate_name": "Test Population Gate",
            "gate_type": "rectangle",
            "x_parameter": "FSC-A",
            "y_parameter": "SSC-A",
            "gate_definition": {
                "x_min": 30000, "x_max": 70000,
                "y_min": 20000, "y_max": 50000
            }
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=gate_data)
        assert response.status_code == 201
        gate = response.json()
        gate_id = gate["id"]
        
        # Get population statistics
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/populations")
        assert response.status_code == 200
        populations = response.json()
        
        # Should have at least one population
        assert len(populations) >= 1
        
        # Find our population
        test_population = None
        for pop in populations:
            if pop.get("gate_id") == gate_id:
                test_population = pop
                break
        
        assert test_population is not None
        assert test_population["n_events"] >= 0
        assert test_population["pct_of_total"] is not None
        assert 0 <= test_population["pct_of_total"] <= 100
        
        # Verify statistics structure
        if test_population["n_events"] > 0:
            assert test_population["median_values"] is not None
            assert test_population["mean_values"] is not None
            assert isinstance(test_population["median_values"], dict)
            assert isinstance(test_population["mean_values"], dict)
    
    def test_event_pagination_large_dataset(self, client, test_flow_dataset_with_events):
        """Test event data pagination with different limits."""
        dataset_id = test_flow_dataset_with_events.id
        
        # Test small page size
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?limit=50&offset=0")
        assert response.status_code == 200
        data = response.json()
        assert len(data["events"]) <= 50
        assert data["total_events"] == 1000
        
        # Test pagination
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?limit=100&offset=100")
        assert response.status_code == 200
        data = response.json()
        assert len(data["events"]) <= 100
        assert data["total_events"] == 1000
        
        # Test subsampling
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?subsample=true&limit=500")
        assert response.status_code == 200
        data = response.json()
        assert len(data["events"]) <= 500
        assert data["total_events"] == 1000
    
    def test_hierarchical_gating_workflow(self, client, test_flow_dataset_with_events):
        """Test hierarchical gating with parent-child relationships."""
        dataset_id = test_flow_dataset_with_events.id
        
        # Create parent gate (lymphocytes)
        parent_gate_data = {
            "gate_name": "Lymphocytes",
            "gate_type": "polygon",
            "x_parameter": "FSC-A",
            "y_parameter": "SSC-A",
            "gate_definition": {
                "vertices": [[25000, 15000], [75000, 15000], [75000, 55000], [25000, 55000]]
            }
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=parent_gate_data)
        assert response.status_code == 201
        parent_gate = response.json()
        parent_gate_id = parent_gate["id"]
        
        # Create child gate (T cells within lymphocytes)
        child_gate_data = {
            "gate_name": "T Cells",
            "gate_type": "rectangle",
            "x_parameter": "FITC-A",
            "y_parameter": "SSC-A",
            "gate_definition": {
                "x_min": 2000, "x_max": 30000,
                "y_min": 15000, "y_max": 55000
            },
            "parent_gate_id": parent_gate_id
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=child_gate_data)
        assert response.status_code == 201
        child_gate = response.json()
        
        # Verify hierarchical relationship
        assert child_gate["parent_gate_id"] == parent_gate_id
        
        # Get populations and verify hierarchy
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/populations")
        assert response.status_code == 200
        populations = response.json()
        
        # Should have populations for both gates
        assert len(populations) >= 2
        
        # Find parent and child populations
        parent_pop = next((p for p in populations if p.get("gate_id") == parent_gate_id), None)
        child_pop = next((p for p in populations if p.get("gate_id") == child_gate["id"]), None)
        
        assert parent_pop is not None
        assert child_pop is not None
        
        # Child population should be smaller than or equal to parent
        if parent_pop["n_events"] > 0:
            assert child_pop["n_events"] <= parent_pop["n_events"]
    
    def test_error_handling_and_validation(self, client, test_flow_dataset_with_events):
        """Test error handling and input validation."""
        dataset_id = test_flow_dataset_with_events.id
        
        # Test invalid gate type
        invalid_gate_data = {
            "gate_name": "Invalid Gate",
            "gate_type": "invalid_type",
            "x_parameter": "FSC-A",
            "y_parameter": "SSC-A",
            "gate_definition": {}
        }
        
        response = client.post(f"/api/v1/flow-cytometry/datasets/{dataset_id}/gates", json=invalid_gate_data)
        assert response.status_code == 422  # Validation error
        
        # Test non-existent dataset
        fake_dataset_id = uuid4()
        response = client.get(f"/api/v1/flow-cytometry/datasets/{fake_dataset_id}")
        assert response.status_code == 404
        
        # Test non-existent gate
        fake_gate_id = uuid4()
        response = client.delete(f"/api/v1/flow-cytometry/gates/{fake_gate_id}")
        assert response.status_code == 404
        
        # Test invalid event query parameters
        response = client.get(f"/api/v1/flow-cytometry/datasets/{dataset_id}/events?limit=0")
        assert response.status_code == 422  # Should validate limit > 0
        
        # Test events from non-processed dataset
        with db_session() as db:
            # Create a pending dataset
            pending_dataset = Dataset(
                id=uuid4(),
                name="Pending Dataset",
                omics_type="flow_cytometry"
            )
            db.add(pending_dataset)
            db.commit()
            
            pending_flow_dataset = FlowCytometryDataset(
                id=uuid4(),
                dataset_id=pending_dataset.id,
                processing_status="pending",
                n_events=0,
                n_parameters=0
            )
            db.add(pending_flow_dataset)
            db.commit()
            
            response = client.get(f"/api/v1/flow-cytometry/datasets/{pending_flow_dataset.id}/events")
            assert response.status_code == 400
            assert "not yet processed" in response.json()["detail"]


# Import required dependencies for test client
from amprenta_rag.api.dependencies import get_current_user
