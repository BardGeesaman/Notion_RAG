"""Integration tests for retrosynthesis API endpoints."""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock
from uuid import uuid4

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


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


class TestRetrosynthesisAPI:
    """Integration tests using real service with mock backend."""

    def test_analyze_valid_smiles(self, client):
        """POST /analyze returns synthesis tree."""
        response = client.post(
            "/api/v1/retrosynthesis/analyze",
            json={"smiles": "CC(=O)Nc1ccccc1", "max_depth": 3}
        )
        assert response.status_code == 200
        data = response.json()
        assert "analysis_id" in data
        assert "tree" in data
        assert data["tree"]["target"] == "CC(=O)Nc1ccccc1"
        assert len(data["tree"]["routes"]) > 0

    def test_analyze_empty_smiles_returns_400(self, client):
        """POST /analyze with empty SMILES returns 400."""
        response = client.post(
            "/api/v1/retrosynthesis/analyze",
            json={"smiles": "", "max_depth": 3}
        )
        assert response.status_code == 400

    def test_analyze_whitespace_smiles_returns_400(self, client):
        """POST /analyze with whitespace SMILES returns 400."""
        response = client.post(
            "/api/v1/retrosynthesis/analyze",
            json={"smiles": "   ", "max_depth": 3}
        )
        assert response.status_code == 400

    def test_get_routes_after_analyze(self, client):
        """GET /routes/{id} returns cached analysis."""
        # First analyze
        analyze_resp = client.post(
            "/api/v1/retrosynthesis/analyze",
            json={"smiles": "CCO", "max_depth": 2}
        )
        assert analyze_resp.status_code == 200
        analysis_id = analyze_resp.json()["analysis_id"]
        
        # Then retrieve
        get_resp = client.get(f"/api/v1/retrosynthesis/routes/{analysis_id}")
        assert get_resp.status_code == 200
        assert get_resp.json()["analysis_id"] == analysis_id

    def test_get_routes_not_found(self, client):
        """GET /routes/{id} returns 404 for unknown ID."""
        response = client.get("/api/v1/retrosynthesis/routes/nonexistent-id")
        assert response.status_code == 404

    def test_building_blocks(self, client):
        """GET /building-blocks returns availability list."""
        response = client.get(
            "/api/v1/retrosynthesis/building-blocks",
            params={"smiles": ["CCO", "CC(=O)O"]}
        )
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 2

    def test_analyze_requires_auth(self):
        """POST /analyze without auth returns 401."""
        # Use client without auth override
        client = TestClient(app)
        response = client.post(
            "/api/v1/retrosynthesis/analyze",
            json={"smiles": "CC", "max_depth": 3}
        )
        assert response.status_code == 401