"""
Unit tests for reports API endpoints.

Tests report generation with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestGenerateReport:
    """Tests for POST /api/v1/reports/generate endpoint."""

    @patch("amprenta_rag.api.routers.reports.generate_report")
    def test_generate_report_success(self, mock_generate):
        """Test successful report generation."""
        entity_id = uuid4()
        
        # Mock report generation
        mock_generate.return_value = "/data/reports/report_123.html"
        
        response = client.post(
            "/api/v1/reports/generate",
            json={
                "entity_type": "program",
                "entity_id": str(entity_id),
                "format": "html",
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["file_path"] == "/data/reports/report_123.html"
        assert "download_url" in data

    @patch("amprenta_rag.api.routers.reports.generate_report")
    def test_generate_report_template_not_found(self, mock_generate):
        """Test report generation with missing template."""
        entity_id = uuid4()
        
        mock_generate.side_effect = FileNotFoundError("Template not found")
        
        response = client.post(
            "/api/v1/reports/generate",
            json={
                "entity_type": "program",
                "entity_id": str(entity_id),
                "format": "html",
            },
        )
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.reports.generate_report")
    def test_generate_report_invalid_entity(self, mock_generate):
        """Test report generation with invalid entity."""
        entity_id = uuid4()
        
        mock_generate.side_effect = ValueError("Invalid entity type")
        
        response = client.post(
            "/api/v1/reports/generate",
            json={
                "entity_type": "invalid_type",
                "entity_id": str(entity_id),
                "format": "html",
            },
        )
        
        assert response.status_code == 400
        assert "invalid" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.reports.generate_report")
    def test_generate_report_generic_error(self, mock_generate):
        """Test report generation with generic error."""
        entity_id = uuid4()
        
        mock_generate.side_effect = RuntimeError("Unexpected error")
        
        response = client.post(
            "/api/v1/reports/generate",
            json={
                "entity_type": "program",
                "entity_id": str(entity_id),
                "format": "html",
            },
        )
        
        assert response.status_code == 500
        assert "failed" in response.json()["detail"].lower()

