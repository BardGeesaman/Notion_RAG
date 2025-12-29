"""
Unit tests for publication extraction API endpoints.

Tests extraction and supplementary file endpoints.
"""

from __future__ import annotations

from unittest.mock import MagicMock
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestExtractExperiments:
    """Tests for POST /api/v1/papers/{paper_id}/extract endpoint."""

    def test_extract_experiments_paper_not_found(self):
        """Test extraction for non-existent paper."""
        paper_id = uuid4()
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(f"/api/v1/papers/{paper_id}/extract")
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()

    def test_extract_experiments_placeholder(self):
        """Test extraction returns placeholder (PDF storage not implemented)."""
        paper_id = uuid4()
        
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_literature
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(f"/api/v1/papers/{paper_id}/extract")
            
            assert response.status_code == 200
            data = response.json()
            assert "experiments_extracted" in data
            assert "extraction_confidence" in data
        finally:
            app.dependency_overrides.clear()


class TestGetPaperExperiments:
    """Tests for GET /api/v1/papers/{paper_id}/experiments endpoint."""

    def test_get_experiments_success(self):
        """Test successful experiment retrieval."""
        paper_id = uuid4()
        exp_id = uuid4()
        
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_exp = MagicMock()
        mock_exp.id = exp_id
        mock_exp.experiment_type = "RNA-seq"
        mock_exp.cell_line = "HeLa"
        mock_exp.treatment = "Drug A"
        mock_exp.concentration = "10 μM"
        mock_exp.timepoint = "24h"
        mock_exp.replicate_count = 3
        mock_exp.measured_entities = ["gene1", "gene2"]
        mock_exp.key_findings = "Significant upregulation"
        mock_exp.extraction_confidence = 85
        
        mock_session = MagicMock()
        
        # First query: verify paper exists
        # Second query: get experiments
        call_count = [0]
        def mock_query_side_effect(*args):
            call_count[0] += 1
            mock_query = MagicMock()
            mock_query.filter.return_value = mock_query
            
            if call_count[0] == 1:
                mock_query.first.return_value = mock_literature
            else:
                mock_query.all.return_value = [mock_exp]
            
            return mock_query
        
        mock_session.query.side_effect = mock_query_side_effect
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/papers/{paper_id}/experiments")
            
            assert response.status_code == 200
            data = response.json()
            assert len(data) == 1
            assert data[0]["experiment_type"] == "RNA-seq"
            assert data[0]["cell_line"] == "HeLa"
        finally:
            app.dependency_overrides.clear()

    def test_get_experiments_paper_not_found(self):
        """Test experiments retrieval for non-existent paper."""
        paper_id = uuid4()
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.get(f"/api/v1/papers/{paper_id}/experiments")
            
            assert response.status_code == 404
        finally:
            app.dependency_overrides.clear()


class TestUploadSupplementary:
    """Tests for POST /api/v1/papers/{paper_id}/supplementary endpoint."""

    def test_upload_supplementary_paper_not_found(self):
        """Test supplementary upload for non-existent paper."""
        paper_id = uuid4()
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            # Create a fake file upload
            files = {"file": ("test.csv", b"Gene,log2FC\nTP53,2.5", "text/csv")}
            response = client.post(
                f"/api/v1/papers/{paper_id}/supplementary",
                files=files,
            )
            
            assert response.status_code == 404
            assert "not found" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()


class TestLinkSupplementaryToDataset:
    """Tests for POST /api/v1/papers/{paper_id}/link-dataset endpoint."""

    def test_link_dataset_success(self):
        """Test successful dataset linking."""
        paper_id = uuid4()
        supp_file_id = uuid4()
        dataset_id = uuid4()
        
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_supp_file = MagicMock()
        mock_supp_file.id = supp_file_id
        mock_supp_file.literature_id = paper_id
        mock_supp_file.linked_dataset_id = None
        
        mock_session = MagicMock()
        
        # Two queries: literature and supplementary file
        call_count = [0]
        def mock_query_side_effect(*args):
            call_count[0] += 1
            mock_query = MagicMock()
            mock_query.filter.return_value = mock_query
            
            if call_count[0] == 1:
                mock_query.first.return_value = mock_literature
            else:
                mock_query.first.return_value = mock_supp_file
            
            return mock_query
        
        mock_session.query.side_effect = mock_query_side_effect
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                f"/api/v1/papers/{paper_id}/link-dataset",
                json={
                    "supplementary_file_id": str(supp_file_id),
                    "dataset_id": str(dataset_id),
                },
            )
            
            assert response.status_code == 200
            data = response.json()
            assert data["linked"] is True
            assert data["supplementary_file_id"] == str(supp_file_id)
            assert data["dataset_id"] == str(dataset_id)
        finally:
            app.dependency_overrides.clear()

    def test_link_dataset_paper_not_found(self):
        """Test linking for non-existent paper."""
        paper_id = uuid4()
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                f"/api/v1/papers/{paper_id}/link-dataset",
                json={
                    "supplementary_file_id": str(uuid4()),
                    "dataset_id": str(uuid4()),
                },
            )
            
            assert response.status_code == 404
        finally:
            app.dependency_overrides.clear()

    def test_link_dataset_supplementary_not_found(self):
        """Test linking when supplementary file not found."""
        paper_id = uuid4()
        
        mock_literature = MagicMock()
        mock_literature.id = paper_id
        
        mock_session = MagicMock()
        
        call_count = [0]
        def mock_query_side_effect(*args):
            call_count[0] += 1
            mock_query = MagicMock()
            mock_query.filter.return_value = mock_query
            
            if call_count[0] == 1:
                mock_query.first.return_value = mock_literature
            else:
                mock_query.first.return_value = None  # Supp file not found
            
            return mock_query
        
        mock_session.query.side_effect = mock_query_side_effect
        
        def mock_get_db():
            yield mock_session
        
        from amprenta_rag.api.dependencies import get_database_session
        
        app.dependency_overrides[get_database_session] = mock_get_db
        
        try:
            response = client.post(
                f"/api/v1/papers/{paper_id}/link-dataset",
                json={
                    "supplementary_file_id": str(uuid4()),
                    "dataset_id": str(uuid4()),
                },
            )
            
            assert response.status_code == 404
            assert "supplementary" in response.json()["detail"].lower()
        finally:
            app.dependency_overrides.clear()

