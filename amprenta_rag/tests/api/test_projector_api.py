"""
Unit tests for projector API endpoints.

Tests projection computation and export.
"""

from __future__ import annotations

from unittest.mock import MagicMock
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestComputeProjection:
    """Tests for POST /api/v1/projector/compute endpoint."""

    def test_compute_pca_with_matrix(self):
        """Test PCA projection with raw matrix."""
        matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        
        response = client.post(
            "/api/v1/projector/compute",
            json={
                "matrix": matrix,
                "algorithm": "pca",
                "n_components": 2,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["algorithm_used"] == "pca"
        assert data["n_samples"] == 3
        assert len(data["coordinates"]) == 3
        assert len(data["coordinates"][0]) == 2

    def test_compute_with_dataset_id(self):
        """Test projection with dataset_id."""
        dataset_id = uuid4()
        
        mock_dataset = MagicMock()
        mock_dataset.id = dataset_id
        mock_dataset.name = "Test Dataset"
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_dataset
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                "/api/v1/projector/compute",
                json={
                    "dataset_id": str(dataset_id),
                    "algorithm": "pca",
                    "n_components": 2,
                },
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["algorithm_used"] == "pca"
        finally:
            app.dependency_overrides.clear()

    def test_compute_missing_input(self):
        """Test compute without dataset_id or matrix."""
        response = client.post(
            "/api/v1/projector/compute",
            json={
                "algorithm": "pca",
                "n_components": 2,
            },
        )
        
        assert response.status_code == 400
        assert "dataset_id or matrix" in response.json()["detail"].lower()

    def test_compute_dataset_not_found(self):
        """Test compute with non-existent dataset."""
        dataset_id = uuid4()
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                "/api/v1/projector/compute",
                json={
                    "dataset_id": str(dataset_id),
                    "algorithm": "pca",
                },
            )
            
            assert response.status_code == 404
        finally:
            app.dependency_overrides.clear()

    def test_compute_unsupported_algorithm(self):
        """Test compute with unsupported algorithm."""
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        
        response = client.post(
            "/api/v1/projector/compute",
            json={
                "matrix": matrix,
                "algorithm": "invalid_algo",
            },
        )
        
        # Should get error (could be 400 or 500 depending on validation point)
        assert response.status_code in [400, 500]
        assert "unsupported" in response.json()["detail"].lower() or "failed" in response.json()["detail"].lower()


class TestListDatasets:
    """Tests for GET /api/v1/projector/datasets endpoint."""

    def test_list_datasets_success(self):
        """Test successful dataset listing."""
        dataset1_id = uuid4()
        dataset2_id = uuid4()
        
        mock_ds1 = MagicMock()
        mock_ds1.id = dataset1_id
        mock_ds1.name = "Dataset 1"
        mock_ds1.omics_type = "transcriptomics"
        mock_ds1.created_at = None
        
        mock_ds2 = MagicMock()
        mock_ds2.id = dataset2_id
        mock_ds2.name = "Dataset 2"
        mock_ds2.omics_type = "proteomics"
        mock_ds2.created_at = None
        
        mock_session = MagicMock()
        mock_session.query.return_value.order_by.return_value.limit.return_value.all.return_value = [
            mock_ds1,
            mock_ds2,
        ]
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get("/api/v1/projector/datasets")
            
            assert response.status_code == 200
            data = response.json()
            assert data["total"] == 2
            assert len(data["datasets"]) == 2
        finally:
            app.dependency_overrides.clear()


class TestExportProjection:
    """Tests for POST /api/v1/projector/export endpoint."""

    def test_export_2d_projection(self):
        """Test exporting 2D projection coordinates."""
        response = client.post(
            "/api/v1/projector/export",
            json={
                "coordinates": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "csv_data" in data
        assert "Sample,X,Y" in data["csv_data"]
        assert "projection_2d.csv" in data["filename"]

    def test_export_3d_projection(self):
        """Test exporting 3D projection coordinates."""
        response = client.post(
            "/api/v1/projector/export",
            json={
                "coordinates": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert "Sample,X,Y,Z" in data["csv_data"]
        assert "projection_3d.csv" in data["filename"]

    def test_export_empty_coordinates(self):
        """Test export with empty coordinates."""
        response = client.post(
            "/api/v1/projector/export",
            json={"coordinates": []},
        )
        
        # Empty coordinates should be handled
        assert response.status_code in [400, 422]  # Either validation or business logic error

