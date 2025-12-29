"""
Unit tests for AutoML API endpoints.

Tests AutoML template listing and launching with mocked dependencies.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestAutoMLTemplatesList:
    """Tests for GET /api/automl/templates endpoint."""

    @patch("amprenta_rag.api.routers.automl._list_templates")
    def test_list_templates_success(self, mock_list_templates):
        """Test successful template list retrieval."""
        mock_list_templates.return_value = [
            {
                "id": "random_forest_classifier",
                "title": "Random Forest Classifier",
                "path": "templates/automl/random_forest_classifier.ipynb",
                "jupyter_url": "http://localhost:8000/notebooks/templates/automl/random_forest_classifier.ipynb",
            },
            {
                "id": "xgboost_regressor",
                "title": "XGBoost Regressor",
                "path": "templates/automl/xgboost_regressor.ipynb",
                "jupyter_url": "http://localhost:8000/notebooks/templates/automl/xgboost_regressor.ipynb",
            },
        ]
        
        response = client.get("/api/automl/templates")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        assert data[0]["id"] == "random_forest_classifier"
        assert data[1]["id"] == "xgboost_regressor"

    @patch("amprenta_rag.api.routers.automl._list_templates")
    def test_list_templates_empty(self, mock_list_templates):
        """Test listing templates when none exist."""
        mock_list_templates.return_value = []
        
        response = client.get("/api/automl/templates")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []


class TestAutoMLLaunch:
    """Tests for POST /api/automl/launch endpoint."""

    @patch("amprenta_rag.api.routers.automl._list_templates")
    def test_launch_template_jupyter_mode(self, mock_list_templates):
        """Test launching template in Jupyter mode."""
        mock_list_templates.return_value = [
            {
                "id": "test_template",
                "title": "Test Template",
                "path": "templates/automl/test_template.ipynb",
                "jupyter_url": "http://localhost:8000/notebooks/templates/automl/test_template.ipynb",
            }
        ]
        
        response = client.post(
            "/api/automl/launch",
            json={
                "template_id": "test_template",
                "params": {"dataset_id": "test123"},
                "run_mode": "jupyter",
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["mode"] == "jupyter"
        assert data["template"]["id"] == "test_template"
        assert data["run_id"] is None

    @patch("amprenta_rag.api.routers.automl._list_templates")
    def test_launch_template_not_found(self, mock_list_templates):
        """Test launching non-existent template."""
        mock_list_templates.return_value = []
        
        response = client.post(
            "/api/automl/launch",
            json={
                "template_id": "nonexistent_template",
                "params": {},
                "run_mode": "jupyter",
            },
        )
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    @patch("amprenta_rag.api.routers.automl._list_templates")
    def test_launch_template_invalid_mode(self, mock_list_templates):
        """Test launching with invalid run mode."""
        mock_list_templates.return_value = [
            {
                "id": "test_template",
                "title": "Test Template",
                "path": "templates/automl/test_template.ipynb",
                "jupyter_url": "http://localhost:8000/notebooks/templates/automl/test_template.ipynb",
            }
        ]
        
        response = client.post(
            "/api/automl/launch",
            json={
                "template_id": "test_template",
                "params": {},
                "run_mode": "invalid_mode",
            },
        )
        
        assert response.status_code == 400
        assert "invalid" in response.json()["detail"].lower()


class TestAutoMLRunsGet:
    """Tests for GET /api/automl/runs/{run_id} endpoint."""

    @patch("amprenta_rag.api.routers.automl._RUNS", {"test-run-id": {
        "id": "test-run-id",
        "template_id": "test_template",
        "status": "success",
        "started_at": 1000.0,
        "ended_at": 1100.0,
        "params": {"dataset_id": "test123"},
        "output_ipynb": "/path/to/output.ipynb",
        "error": None,
    }})
    def test_get_run_success(self):
        """Test successful run retrieval."""
        response = client.get("/api/automl/runs/test-run-id")
        
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == "test-run-id"
        assert data["status"] == "success"
        assert data["template_id"] == "test_template"

    @patch("amprenta_rag.api.routers.automl._RUNS", {})
    def test_get_run_not_found(self):
        """Test getting non-existent run."""
        response = client.get("/api/automl/runs/nonexistent-run-id")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestAutoMLRunsList:
    """Tests for GET /api/automl/runs endpoint."""

    @patch("amprenta_rag.api.routers.automl._RUNS", {
        "run1": {"id": "run1", "status": "success"},
        "run2": {"id": "run2", "status": "running"},
        "run3": {"id": "run3", "status": "failed"},
    })
    def test_list_runs_success(self):
        """Test successful runs list retrieval."""
        response = client.get("/api/automl/runs")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 3
        assert any(r["id"] == "run1" for r in data)
        assert any(r["id"] == "run2" for r in data)
        assert any(r["id"] == "run3" for r in data)

    @patch("amprenta_rag.api.routers.automl._RUNS", {})
    def test_list_runs_empty(self):
        """Test listing runs when none exist."""
        response = client.get("/api/automl/runs")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []

