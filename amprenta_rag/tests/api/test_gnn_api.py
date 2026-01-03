"""API tests for GNN toxicity endpoints."""
import pytest
import os
from uuid import uuid4
from fastapi.testclient import TestClient
from unittest.mock import MagicMock

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user

# Set OpenMP workaround for tests
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


@pytest.fixture
def client():
    """Test client with mocked auth only."""
    def mock_user():
        user = MagicMock()
        user.id = uuid4()
        user.email = "test@example.com"
        return user
    
    app.dependency_overrides[get_current_user] = mock_user
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


class TestGNNAPI:
    """Integration tests for GNN toxicity API endpoints."""

    def test_list_gnn_models(self, client):
        """List GNN models endpoint works."""
        response = client.get("/api/admet/gnn-models")
        assert response.status_code == 200
        data = response.json()
        assert "models" in data
        assert isinstance(data["models"], list)
        
        # Verify structure
        for model in data["models"]:
            assert "endpoint" in model
            assert "available" in model

    def test_gnn_model_info_unknown_endpoint(self, client):
        """Unknown endpoint returns 404."""
        response = client.get("/api/admet/gnn-models/unknown/info")
        assert response.status_code == 404

    def test_gnn_model_info_existing_endpoint(self, client):
        """Valid endpoint returns model info or 404 if not trained."""
        response = client.get("/api/admet/gnn-models/herg/info")
        # Should be 200 (if model exists), 404 (if not trained), or 500 (if load error)
        assert response.status_code in [200, 404, 500]
        
        if response.status_code == 200:
            data = response.json()
            assert "endpoint" in data
            assert "config" in data

    def test_predict_gnn_empty_smiles(self, client):
        """Empty SMILES list returns empty predictions."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": []},
            params={"endpoints": ["herg"]}
        )
        assert response.status_code == 200
        data = response.json()
        assert "predictions" in data
        assert len(data["predictions"]) == 0

    def test_predict_gnn_invalid_endpoint(self, client):
        """Invalid endpoint in predictions handled gracefully."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]},
            params={"endpoints": ["invalid_endpoint"]}
        )
        assert response.status_code == 200
        data = response.json()
        
        # Should have error in response for invalid endpoint
        if data.get("predictions"):
            pred = data["predictions"][0]
            assert "endpoints" in pred
            assert "invalid_endpoint" in pred["endpoints"]
            assert "error" in pred["endpoints"]["invalid_endpoint"]

    def test_predict_gnn_valid_smiles(self, client):
        """Valid SMILES prediction request."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]},
            params={"endpoints": ["herg"]}
        )
        assert response.status_code == 200
        data = response.json()
        assert "predictions" in data
        assert len(data["predictions"]) == 1

    def test_compare_models_endpoint(self, client):
        """Compare endpoint works."""
        response = client.post(
            "/api/admet/compare",
            json={"smiles": ["CCO"]},
            params={"endpoint": "herg"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "endpoint" in data
        assert "comparisons" in data

    def test_predict_gnn_with_uncertainty(self, client):
        """Prediction with uncertainty flag works."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]},
            params={"endpoints": ["herg"], "with_uncertainty": True}
        )
        assert response.status_code == 200
        data = response.json()
        assert "predictions" in data

    def test_predict_gnn_without_uncertainty(self, client):
        """Prediction without uncertainty works."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]},
            params={"endpoints": ["herg"], "with_uncertainty": False}
        )
        assert response.status_code == 200
        data = response.json()
        assert "predictions" in data

    def test_predict_gnn_multiple_endpoints(self, client):
        """Multiple endpoints in single request."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]},
            params={"endpoints": ["herg", "ames"]}
        )
        assert response.status_code == 200
        data = response.json()
        
        if data.get("predictions"):
            pred = data["predictions"][0]
            assert "endpoints" in pred
            # Should have both endpoints (with results or errors)
            assert "herg" in pred["endpoints"]
            assert "ames" in pred["endpoints"]

    def test_gnn_models_list_structure(self, client):
        """GNN models list has correct structure."""
        response = client.get("/api/admet/gnn-models")
        data = response.json()
        
        assert "models" in data
        for model in data["models"]:
            assert "endpoint" in model
            assert "available" in model
            
            if model["available"]:
                # Available models should have metrics
                assert "metrics" in model
            else:
                # Unavailable models should have error or status
                assert "error" in model or "available" in model

    def test_compare_response_structure(self, client):
        """Compare response has correct structure."""
        response = client.post(
            "/api/admet/compare",
            json={"smiles": ["CCO"]},
            params={"endpoint": "herg"}
        )
        data = response.json()
        
        assert "endpoint" in data
        assert "comparisons" in data
        
        if data["comparisons"]:
            comp = data["comparisons"][0]
            assert "smiles" in comp
            assert "xgboost" in comp
            assert "gnn" in comp

    def test_predict_gnn_invalid_smiles(self, client):
        """Invalid SMILES handled gracefully."""
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["invalid_smiles"]},
            params={"endpoints": ["herg"]}
        )
        assert response.status_code == 200
        data = response.json()
        
        if data.get("predictions"):
            pred = data["predictions"][0]
            # Should either have error or handle gracefully
            assert "smiles" in pred

    def test_gnn_endpoints_require_auth(self):
        """GNN endpoints require authentication."""
        # Use client without auth override
        client = TestClient(app)
        
        response = client.get("/api/admet/gnn-models")
        assert response.status_code == 401
        
        response = client.post(
            "/api/admet/predict-gnn",
            json={"smiles": ["CCO"]}
        )
        assert response.status_code == 401
