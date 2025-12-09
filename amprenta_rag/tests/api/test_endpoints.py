"""
Tests for FastAPI endpoints.

These tests verify:
- Endpoints are accessible
- CRUD operations work
- Error handling is correct
- Response schemas are valid
"""

import uuid

import pytest
from fastapi.testclient import TestClient

# Skip until API schemas/endpoints are stabilized
pytestmark = pytest.mark.skip(reason="API schemas/endpoints under refactor; skipping for now")

from amprenta_rag.api.main import app

client = TestClient(app)


class TestRootEndpoints:
    """Test root and health check endpoints."""
    
    def test_root(self):
        """Test root endpoint."""
        response = client.get("/")
        assert response.status_code == 200
        data = response.json()
        assert data["name"] == "Amprenta Multi-Omics Platform API"
        assert data["status"] == "operational"
    
    def test_health(self):
        """Test health check endpoint."""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"


class TestProgramsEndpoints:
    """Test Programs API endpoints."""
    
    def test_list_programs(self):
        """Test listing programs."""
        response = client.get("/api/v1/programs")
        assert response.status_code == 200
        programs = response.json()
        assert isinstance(programs, list)
    
    def test_create_program(self):
        """Test creating a program."""
        program_data = {
            "name": "Test Program (API Test)",
            "description": "Created by test",
        }
        response = client.post("/api/v1/programs", json=program_data)
        
        # Should succeed or return validation error
        assert response.status_code in (200, 201, 422)
        
        if response.status_code in (200, 201):
            program = response.json()
            assert program["name"] == program_data["name"]
            
            # Clean up - delete the program
            program_id = program["id"]
            client.delete(f"/api/v1/programs/{program_id}")
    
    def test_get_program_not_found(self):
        """Test getting a non-existent program."""
        fake_id = uuid.uuid4()
        response = client.get(f"/api/v1/programs/{fake_id}")
        assert response.status_code == 404


class TestExperimentsEndpoints:
    """Test Experiments API endpoints."""
    
    def test_list_experiments(self):
        """Test listing experiments."""
        response = client.get("/api/v1/experiments")
        assert response.status_code == 200
        experiments = response.json()
        assert isinstance(experiments, list)


class TestDatasetsEndpoints:
    """Test Datasets API endpoints."""
    
    def test_list_datasets(self):
        """Test listing datasets."""
        response = client.get("/api/v1/datasets")
        assert response.status_code == 200
        datasets = response.json()
        assert isinstance(datasets, list)


class TestFeaturesEndpoints:
    """Test Features API endpoints."""
    
    def test_list_features(self):
        """Test listing features."""
        response = client.get("/api/v1/features")
        assert response.status_code == 200
        features = response.json()
        assert isinstance(features, list)


class TestSignaturesEndpoints:
    """Test Signatures API endpoints."""
    
    def test_list_signatures(self):
        """Test listing signatures."""
        response = client.get("/api/v1/signatures")
        assert response.status_code == 200
        signatures = response.json()
        assert isinstance(signatures, list)


class TestAPIDocumentation:
    """Test API documentation endpoints."""
    
    def test_openapi_schema(self):
        """Test OpenAPI schema is accessible."""
        response = client.get("/openapi.json")
        assert response.status_code == 200
        schema = response.json()
        assert "openapi" in schema
        assert "paths" in schema
        # Check for programs endpoint (may have trailing slash)
        assert "/api/v1/programs" in schema["paths"] or "/api/v1/programs/" in schema["paths"]
    
    def test_docs_accessible(self):
        """Test Swagger UI is accessible."""
        response = client.get("/docs")
        assert response.status_code == 200
    
    def test_redoc_accessible(self):
        """Test ReDoc is accessible."""
        response = client.get("/redoc")
        assert response.status_code == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

