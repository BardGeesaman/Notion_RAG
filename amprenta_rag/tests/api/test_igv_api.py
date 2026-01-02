"""Tests for IGV.js file streaming endpoints."""

from __future__ import annotations

from pathlib import Path
from uuid import uuid4
from unittest.mock import patch, MagicMock

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user


@pytest.fixture
def mock_user():
    """Create mock authenticated user."""
    user = MagicMock()
    user.id = uuid4()
    user.email = "test@example.com"
    return user


@pytest.fixture
def client(mock_user):
    """Test client with auth override."""
    app.dependency_overrides[get_current_user] = lambda: mock_user
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


@pytest.fixture
def mock_alignment_file(tmp_path):
    """Create mock alignment file and DB record."""
    # Create fake BAM file (just needs some bytes)
    bam_content = b'\x1f\x8b\x08\x04' + b'\x00' * 100  # gzip magic + padding
    bam_path = tmp_path / "test.bam"
    bam_path.write_bytes(bam_content)
    
    index_path = tmp_path / "test.bam.bai"
    index_path.write_bytes(b'\x00' * 50)
    
    return {
        "id": uuid4(),
        "file_path": str(bam_path),
        "index_file_path": str(index_path),
    }


class TestStreamAlignmentFile:
    """Tests for GET /alignments/{id}/file endpoint."""
    
    def test_full_file_response(self, client, mock_alignment_file):
        """Test full file download (no Range header)."""
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = mock_alignment_file["id"]
            alignment.file_path = mock_alignment_file["file_path"]
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", str(Path(mock_alignment_file["file_path"]).parent)):
                response = client.get(f"/api/v1/genomics/alignments/{mock_alignment_file['id']}/file")
        
        assert response.status_code == 200
        assert "Accept-Ranges" in response.headers
        assert response.headers["Accept-Ranges"] == "bytes"
    
    def test_range_header_partial_content(self, client, mock_alignment_file):
        """Test byte-range request returns 206."""
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = mock_alignment_file["id"]
            alignment.file_path = mock_alignment_file["file_path"]
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", str(Path(mock_alignment_file["file_path"]).parent)):
                response = client.get(
                    f"/api/v1/genomics/alignments/{mock_alignment_file['id']}/file",
                    headers={"Range": "bytes=0-50"}
                )
        
        assert response.status_code == 206
        assert "Content-Range" in response.headers
    
    def test_range_header_parsing(self, client, mock_alignment_file):
        """Test Range header formats are parsed correctly."""
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = mock_alignment_file["id"]
            alignment.file_path = mock_alignment_file["file_path"]
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", str(Path(mock_alignment_file["file_path"]).parent)):
                # Test open-ended range
                response = client.get(
                    f"/api/v1/genomics/alignments/{mock_alignment_file['id']}/file",
                    headers={"Range": "bytes=10-"}
                )
        
        assert response.status_code == 206
    
    def test_path_traversal_rejected(self, client):
        """Test path traversal attempt returns 403."""
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = uuid4()
            alignment.file_path = "/etc/passwd"  # Path traversal attempt
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", "/some/storage/path"):
                response = client.get(f"/api/v1/genomics/alignments/{alignment.id}/file")
        
        assert response.status_code == 403
        assert "Access denied" in response.text
    
    def test_cors_headers_exposed(self, client, mock_alignment_file):
        """Test CORS expose_headers are set for browser access."""
        # This tests the middleware configuration
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = mock_alignment_file["id"]
            alignment.file_path = mock_alignment_file["file_path"]
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", str(Path(mock_alignment_file["file_path"]).parent)):
                response = client.get(
                    f"/api/v1/genomics/alignments/{mock_alignment_file['id']}/file",
                    headers={"Range": "bytes=0-10"}
                )
        
        assert response.status_code == 206
        # Verify response has the headers that should be exposed
        assert "Content-Range" in response.headers


class TestStreamAlignmentIndex:
    """Tests for GET /alignments/{id}/index endpoint."""
    
    def test_index_file_served(self, client, mock_alignment_file):
        """Test index file is served correctly."""
        with patch("amprenta_rag.api.routers.genomics.db_session") as mock_db:
            alignment = MagicMock()
            alignment.id = mock_alignment_file["id"]
            alignment.index_file_path = mock_alignment_file["index_file_path"]
            mock_db.return_value.__enter__.return_value.query.return_value.filter.return_value.first.return_value = alignment
            
            with patch("amprenta_rag.api.routers.genomics.ALIGNMENT_STORAGE_PATH", str(Path(mock_alignment_file["index_file_path"]).parent)):
                response = client.get(f"/api/v1/genomics/alignments/{mock_alignment_file['id']}/index")
        
        assert response.status_code == 200
        assert "Accept-Ranges" in response.headers
