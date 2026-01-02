"""Tests for genomics API endpoints."""

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


class TestVCFUpload:
    """Tests for VCF upload endpoint."""

    def test_vcf_upload_rejects_large_file(self, client):
        """Files > 50MB are rejected with 413."""
        # Create fake large content (just check the size validation logic)
        large_content = b"x" * (51 * 1024 * 1024)  # 51MB
        
        response = client.post(
            "/api/v1/genomics/variants/upload",
            files={"file": ("test.vcf", large_content, "text/plain")},
        )
        
        assert response.status_code == 413

    def test_vcf_upload_rejects_invalid_type(self, client):
        """Non-VCF files are rejected with 415."""
        response = client.post(
            "/api/v1/genomics/variants/upload",
            files={"file": ("test.txt", b"not a vcf", "text/plain")},
        )
        
        assert response.status_code == 415

    @patch("amprenta_rag.api.routers.genomics.parse_vcf_bytes")
    def test_vcf_upload_preview_mode(self, mock_parse, client):
        """Preview mode returns variants without saving."""
        mock_parse.return_value = [
            {"chromosome": "chr1", "position": 123, "gene": "TEST"}
        ]
        
        response = client.post(
            "/api/v1/genomics/variants/upload?preview_only=true",
            files={"file": ("test.vcf", b"##fileformat=VCF", "text/plain")},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["variants_count"] == 1
        assert len(data["variants"]) == 1


class TestENASearch:
    """Tests for ENA search endpoint."""

    @patch("amprenta_rag.api.routers.genomics.ENARepository")
    def test_ena_search_returns_results(self, mock_repo_class, client):
        """Search endpoint returns ENA results."""
        mock_repo = MagicMock()
        mock_repo.search_studies.return_value = ["ERR123456"]
        mock_repo.fetch_study_metadata.return_value = MagicMock(
            title="Test Study",
            organism=["Homo sapiens"],
            platform="Illumina",
            raw_metadata={"study_accession": "PRJNA123", "fastq_ftp": "ftp://a.fq;ftp://b.fq"}
        )
        mock_repo_class.return_value = mock_repo
        
        response = client.get("/api/v1/genomics/ena/search?q=cancer")
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["total"] >= 0


class TestENAIngest:
    """Tests for ENA ingest endpoint."""

    def test_ena_ingest_creates_job(self, client):
        """Ingest endpoint accepts study IDs."""
        response = client.post(
            "/api/v1/genomics/ena/ingest",
            json={"study_ids": ["ERR123456", "ERR789012"]},
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "queued" in data["message"].lower() or "2" in data["message"]
