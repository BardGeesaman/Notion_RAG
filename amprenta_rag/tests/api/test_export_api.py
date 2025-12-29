"""Tests for export API endpoints."""

from __future__ import annotations

from unittest.mock import patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestExportDataset:
    """Tests for GET /api/v1/export/dataset/{id} endpoint."""

    @patch("amprenta_rag.api.routers.export.export_dataset")
    def test_export_dataset_csv(self, mock_export):
        """Test dataset CSV export."""
        dataset_id = uuid4()
        mock_export.return_value = b"name,value\ngene1,2.5"
        
        response = client.get(f"/api/v1/export/dataset/{dataset_id}?format=csv")
        
        assert response.status_code == 200
        assert response.headers["content-type"] == "text/csv; charset=utf-8"


class TestExportExperiment:
    """Tests for GET /api/v1/export/experiment/{id} endpoint."""

    @patch("amprenta_rag.api.routers.export.export_experiment")
    def test_export_experiment_json(self, mock_export):
        """Test experiment JSON export."""
        experiment_id = uuid4()
        mock_export.return_value = b'{"experiment_id": "test"}'
        
        response = client.get(f"/api/v1/export/experiment/{experiment_id}?format=json")
        
        assert response.status_code == 200
        assert response.headers["content-type"] == "application/json"


class TestExportCompounds:
    """Tests for POST /api/v1/export/compounds endpoint."""

    @patch("amprenta_rag.api.routers.export.export_compounds")
    def test_export_compounds_csv(self, mock_export):
        """Test compounds CSV export."""
        mock_export.return_value = b"compound_id,smiles\nCMPD001,CCO"
        
        response = client.post(
            "/api/v1/export/compounds",
            json={
                "compound_ids": [str(uuid4())],
                "format": "csv",
            },
        )
        
        assert response.status_code == 200
        assert "text/csv" in response.headers["content-type"]


class TestCreatePackage:
    """Tests for POST /api/v1/export/package endpoint."""

    @patch("amprenta_rag.api.routers.export.create_export_package")
    def test_create_package_zip(self, mock_create):
        """Test package creation."""
        mock_create.return_value = b"PK\x03\x04"  # ZIP magic number
        
        response = client.post(
            "/api/v1/export/package",
            json={
                "items": ["dataset1", "dataset2"],
                "format": "zip",
            },
        )
        
        assert response.status_code == 200
        assert response.headers["content-type"] == "application/zip"

