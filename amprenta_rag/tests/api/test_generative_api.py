"""Tests for generative chemistry API endpoints."""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user
from amprenta_rag.api.routers.generative import get_generative_service
from amprenta_rag.ml.generative.service import GenerativeModelNotFoundError


class TestGenerativeAPI:
    """Test generative chemistry API endpoints."""
    
    @pytest.fixture
    def mock_user(self):
        """Mock authenticated user."""
        return {"id": uuid4(), "username": f"testuser_{uuid4().hex[:8]}"}
    
    @pytest.fixture
    def mock_service(self):
        """Mock generative chemistry service."""
        service = MagicMock()
        
        # Default mock responses
        service.get_model_info.return_value = {
            "name": "test_model",
            "version": "1.0.0",
            "latent_dim": 256,
            "vocab_size": 52,
            "status": "loaded",
        }
        
        service.sample.return_value = {
            "molecules": [
                {"smiles": "CCO", "properties": {"mw": 46.07}, "score": None},
                {"smiles": "CCC", "properties": {"mw": 44.10}, "score": None},
            ],
            "count": 2
        }
        
        service.interpolate.return_value = [
            {"smiles": "CCO", "properties": {}, "score": None, "step": 0},
            {"smiles": "CC(C)O", "properties": {}, "score": None, "step": 1},
            {"smiles": "CCC", "properties": {}, "score": None, "step": 2},
        ]
        
        service.optimize.return_value = {
            "optimized": [
                {"smiles": "CC(C)O", "properties": {"logp": 0.1}, "score": 0.8, "iteration": 5}
            ],
            "seed_properties": {"logp": 0.0},
            "best_score": 0.8,
            "iterations_completed": 10,
        }
        
        service.scaffold_hop.return_value = {
            "scaffold": "c1ccccc1",
            "analogs": [
                {"smiles": "Cc1ccccc1", "properties": {}, "score": None},
                {"smiles": "CCc1ccccc1", "properties": {}, "score": None},
            ],
            "n_generated": 2,
        }
        
        return service
    
    @pytest.fixture
    def client(self, mock_user, mock_service):
        """Test client with mocked dependencies."""
        app.dependency_overrides[get_current_user] = lambda: mock_user
        app.dependency_overrides[get_generative_service] = lambda: mock_service
        
        try:
            client = TestClient(app)
            yield client
        finally:
            # Cleanup
            app.dependency_overrides.clear()
    
    def test_sample_default(self, client, mock_service):
        """Test default molecule sampling."""
        response = client.post("/api/v1/generative/sample", json={
            "n_samples": 2,
            "temperature": 1.0,
            "max_length": 100,
        })
        
        assert response.status_code == 200
        data = response.json()
        
        assert "molecules" in data
        assert "count" in data
        assert "model_info" in data
        assert data["count"] == 2
        assert len(data["molecules"]) == 2
        
        # Verify service was called with correct parameters
        mock_service.sample.assert_called_once_with(
            n_samples=2,
            temperature=1.0,
            max_length=100,
        )
    
    def test_sample_custom_count(self, client, mock_service):
        """Test sampling with custom molecule count."""
        # Configure mock for different count
        mock_service.sample.return_value = [
            {"smiles": f"C{'C' * i}O", "properties": {}, "score": None}
            for i in range(5)
        ]
        
        response = client.post("/api/v1/generative/sample", json={
            "n_samples": 5,
            "temperature": 0.8,
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["count"] == 5
        
        mock_service.sample.assert_called_once_with(
            n_samples=5,
            temperature=0.8,
            max_length=100,  # Default value
        )
    
    def test_sample_invalid_params(self, client):
        """Test sampling with invalid parameters."""
        # Test invalid n_samples (too high)
        response = client.post("/api/v1/generative/sample", json={
            "n_samples": 200,  # Exceeds limit of 100
        })
        assert response.status_code == 422  # Validation error
        
        # Test invalid temperature (too low)
        response = client.post("/api/v1/generative/sample", json={
            "temperature": 0.05,  # Below minimum of 0.1
        })
        assert response.status_code == 422
    
    def test_interpolate_valid(self, client, mock_service):
        """Test valid molecular interpolation."""
        response = client.post("/api/v1/generative/interpolate", json={
            "smiles_start": "CCO",
            "smiles_end": "CCC",
            "steps": 3,
            "interpolation_type": "linear",
        })
        
        assert response.status_code == 200
        data = response.json()
        
        assert "molecules" in data
        assert "start_smiles" in data
        assert "end_smiles" in data
        assert "steps" in data
        assert data["start_smiles"] == "CCO"
        assert data["end_smiles"] == "CCC"
        assert data["steps"] == 3
        assert len(data["molecules"]) == 3
        
        # Check step indices
        steps = [mol["step"] for mol in data["molecules"]]
        assert steps == [0, 1, 2]
        
        mock_service.interpolate.assert_called_once_with(
            smiles_start="CCO",
            smiles_end="CCC",
            steps=3,
            interpolation_type="linear",
        )
    
    def test_interpolate_invalid_smiles(self, client):
        """Test interpolation with invalid SMILES."""
        # Empty SMILES
        response = client.post("/api/v1/generative/interpolate", json={
            "smiles_start": "",
            "smiles_end": "CCC",
            "steps": 3,
        })
        assert response.status_code == 422  # Validation error
        
        # Invalid steps (too few)
        response = client.post("/api/v1/generative/interpolate", json={
            "smiles_start": "CCO",
            "smiles_end": "CCC",
            "steps": 1,  # Below minimum of 2
        })
        assert response.status_code == 422
    
    def test_optimize_basic(self, client, mock_service):
        """Test basic molecular optimization."""
        response = client.post("/api/v1/generative/optimize", json={
            "seed_smiles": "CCO",
            "constraints": [
                {
                    "name": "logp",
                    "min_value": 0.0,
                    "max_value": 2.0,
                    "weight": 1.0,
                }
            ],
            "n_iterations": 10,
        })
        
        assert response.status_code == 200
        data = response.json()
        
        assert "optimized" in data
        assert "seed_smiles" in data
        assert "seed_properties" in data
        assert "best_score" in data
        assert "iterations_completed" in data
        assert data["seed_smiles"] == "CCO"
        assert data["best_score"] == 0.8
        assert len(data["optimized"]) == 1
        
        # Verify service call
        mock_service.optimize.assert_called_once()
        call_args = mock_service.optimize.call_args[1]
        assert call_args["seed_smiles"] == "CCO"
        assert len(call_args["constraints"]) == 1
        assert call_args["constraints"][0]["name"] == "logp"
    
    def test_optimize_with_constraints(self, client, mock_service):
        """Test optimization with multiple constraints."""
        constraints = [
            {
                "name": "logp",
                "min_value": 0.0,
                "max_value": 3.0,
                "target_value": 1.5,
                "weight": 2.0,
            },
            {
                "name": "mw",
                "min_value": 100.0,
                "max_value": 500.0,
                "weight": 1.0,
            },
        ]
        
        response = client.post("/api/v1/generative/optimize", json={
            "seed_smiles": "CCO",
            "constraints": constraints,
            "n_iterations": 50,
            "n_samples_per_iter": 20,
            "learning_rate": 0.2,
            "temperature": 1.5,
        })
        
        assert response.status_code == 200
        
        # Verify all parameters passed correctly
        call_args = mock_service.optimize.call_args[1]
        assert call_args["n_iterations"] == 50
        assert call_args["n_samples_per_iter"] == 20
        assert call_args["learning_rate"] == 0.2
        assert call_args["temperature"] == 1.5
        assert len(call_args["constraints"]) == 2
    
    def test_scaffold_hop(self, client, mock_service):
        """Test scaffold hopping."""
        response = client.post("/api/v1/generative/scaffold-hop", json={
            "smiles": "CCc1ccccc1",
            "n_analogs": 2,
            "preserve_scaffold": False,
            "similarity_threshold": 0.7,
        })
        
        assert response.status_code == 200
        data = response.json()
        
        assert "scaffold" in data
        assert "analogs" in data
        assert "input_smiles" in data
        assert "n_generated" in data
        assert data["scaffold"] == "c1ccccc1"
        assert data["input_smiles"] == "CCc1ccccc1"
        assert data["n_generated"] == 2
        assert len(data["analogs"]) == 2
        
        mock_service.scaffold_hop.assert_called_once_with(
            smiles="CCc1ccccc1",
            n_analogs=2,
            preserve_scaffold=False,
            similarity_threshold=0.7,
        )
    
    def test_list_models(self, client, mock_service):
        """Test listing available models."""
        response = client.get("/api/v1/generative/models")
        
        assert response.status_code == 200
        data = response.json()
        
        assert "models" in data
        assert "count" in data
        assert "default_model" in data
        assert data["count"] == 1
        assert len(data["models"]) == 1
        
        model = data["models"][0]
        assert model["name"] == "test_model"
        assert model["version"] == "1.0.0"
        assert model["latent_dim"] == 256
        assert model["vocab_size"] == 52
        assert model["status"] == "loaded"
        
        mock_service.get_model_info.assert_called_once()
    
    def test_model_not_found_error(self, client, mock_service):
        """Test handling of model not found errors."""
        # Configure service to raise model not found error
        mock_service.sample.side_effect = GenerativeModelNotFoundError("Model not found")
        
        response = client.post("/api/v1/generative/sample", json={
            "n_samples": 2,
        })
        
        assert response.status_code == 503  # Service unavailable
        data = response.json()
        assert "detail" in data
        assert "model not available" in data["detail"].lower()
    
    def test_invalid_constraints_error(self, client):
        """Test optimization with invalid constraints."""
        # Empty constraints list
        response = client.post("/api/v1/generative/optimize", json={
            "seed_smiles": "CCO",
            "constraints": [],  # Empty list not allowed
        })
        assert response.status_code == 422  # Validation error
        
        # Invalid constraint name
        response = client.post("/api/v1/generative/optimize", json={
            "seed_smiles": "CCO",
            "constraints": [
                {
                    "name": "",  # Empty name not allowed
                    "min_value": 0.0,
                    "max_value": 1.0,
                }
            ],
        })
        assert response.status_code == 422
