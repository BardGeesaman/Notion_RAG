"""Tests for alignment API endpoints."""

from unittest.mock import patch, MagicMock
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


@pytest.fixture
def client():
    """Create test client with auth override."""
    from amprenta_rag.api.dependencies import get_current_user
    
    mock_user = MagicMock()
    mock_user.id = uuid4()
    mock_user.email = "test@test.com"
    
    app.dependency_overrides[get_current_user] = lambda: mock_user
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


class TestAlignmentUpload:
    """Tests for alignment upload endpoint."""

    def test_upload_rejects_invalid_format(self, client):
        """Non-BAM/CRAM files are rejected."""
        response = client.post(
            "/api/v1/genomics/alignments/upload",
            files={"file": ("test.txt", b"not a bam file", "text/plain")},
        )
        
        assert response.status_code == 415

    def test_upload_rejects_large_file(self, client):
        """Files > 500MB are rejected."""
        # Create mock large content header
        large_content = b"x" * (501 * 1024 * 1024)
        
        response = client.post(
            "/api/v1/genomics/alignments/upload",
            files={"file": ("test.bam", large_content, "application/octet-stream")},
        )
        
        assert response.status_code == 413


class TestAlignmentList:
    """Tests for alignment list endpoint."""

    @patch("amprenta_rag.api.routers.genomics.db_session")
    def test_list_alignments_returns_empty(self, mock_db, client):
        """List returns empty array when no alignments."""
        mock_session = MagicMock()
        mock_session.query.return_value.order_by.return_value.limit.return_value.all.return_value = []
        mock_db.return_value.__enter__.return_value = mock_session
        
        response = client.get("/api/v1/genomics/alignments")
        
        assert response.status_code == 200
        assert response.json() == []


class TestAlignmentReads:
    """Tests for alignment reads endpoint."""

    @patch("amprenta_rag.api.routers.genomics.db_session")
    def test_reads_requires_index(self, mock_db, client):
        """Reads endpoint requires index file."""
        mock_alignment = MagicMock()
        mock_alignment.id = uuid4()
        mock_alignment.has_index = False
        mock_alignment.file_path = "/fake/path.bam"
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_alignment
        mock_db.return_value.__enter__.return_value = mock_session
        
        response = client.get(
            f"/api/v1/genomics/alignments/{mock_alignment.id}/reads",
            params={"region": "chr1:1000-2000"},
        )
        
        assert response.status_code == 400
        assert "index" in response.text.lower()
