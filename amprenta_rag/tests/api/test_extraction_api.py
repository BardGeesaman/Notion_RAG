"""API tests for extraction endpoints."""

import pytest
from unittest.mock import patch, MagicMock
from uuid import uuid4

from fastapi.testclient import TestClient
from amprenta_rag.api.main import app


@pytest.fixture
def client():
    """Create test client with auth override."""
    mock_user = MagicMock()
    mock_user.id = uuid4()
    mock_user.email = "test@example.com"
    
    from amprenta_rag.api.dependencies import get_current_user
    app.dependency_overrides[get_current_user] = lambda: mock_user
    
    try:
        yield TestClient(app)
    finally:
        app.dependency_overrides.clear()


class TestOCREndpoint:
    """Tests for POST /extraction/ocr endpoint."""

    def test_ocr_rejects_unsupported_file_type(self, client):
        """Unsupported file types return 415."""
        files = {"file": ("test.txt", b"text content", "text/plain")}
        response = client.post("/api/extraction/ocr", files=files)
        assert response.status_code == 415
        assert "Unsupported file type" in response.json()["detail"]

    def test_ocr_rejects_large_file(self, client):
        """Files over 10MB return 413."""
        # Create a file larger than 10MB
        large_content = b"x" * (10_000_001)
        files = {"file": ("large.png", large_content, "image/png")}
        response = client.post("/api/extraction/ocr", files=files)
        assert response.status_code == 413
        assert "too large" in response.json()["detail"].lower()

    @patch("amprenta_rag.extraction.ocr_service.OCRService")
    def test_ocr_success(self, mock_ocr_class, client):
        """Valid file returns extracted text."""
        mock_ocr = MagicMock()
        mock_ocr.extract_from_image.return_value = "Extracted test text"
        mock_ocr_class.return_value = mock_ocr
        
        files = {"file": ("test.png", b"fake image content", "image/png")}
        response = client.post("/api/extraction/ocr?language=eng", files=files)
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["text"] == "Extracted test text"


class TestScrapeEndpoint:
    """Tests for POST /extraction/scrape endpoint."""

    @patch("amprenta_rag.extraction.web_scraper.WebScraper")
    def test_scrape_returns_content(self, mock_scraper_class, client):
        """Successful scrape returns content."""
        mock_result = MagicMock()
        mock_result.success = True
        mock_result.url = "https://example.com"
        mock_result.title = "Test Title"
        mock_result.author = "Test Author"
        mock_result.content = "Test content"
        mock_result.word_count = 2
        mock_result.error = None
        
        mock_scraper = MagicMock()
        mock_scraper.extract_from_url.return_value = mock_result
        mock_scraper_class.return_value = mock_scraper
        
        response = client.post("/api/extraction/scrape", json={"url": "https://example.com"})
        
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["title"] == "Test Title"


class TestNormalizeEndpoint:
    """Tests for POST /extraction/normalize endpoint."""

    def test_normalize_rejects_invalid_type(self, client):
        """Invalid entity type returns 400."""
        response = client.post("/api/extraction/normalize", json={
            "entity_type": "invalid",
            "name": "test"
        })
        assert response.status_code == 400
        assert "Invalid entity_type" in response.json()["detail"]


class TestSSRFPrevention:
    """SSRF prevention tests for /scrape endpoint."""

    def test_scrape_rejects_private_ip(self, client):
        """Private IPv4 should be blocked."""
        response = client.post(
            "/api/extraction/scrape",
            json={"url": "http://10.0.0.1/"},
        )
        assert response.status_code == 400
        assert "blocked" in response.json()["detail"].lower()

    def test_scrape_rejects_localhost(self, client):
        """Localhost should be blocked."""
        response = client.post(
            "/api/extraction/scrape",
            json={"url": "http://localhost/"},
        )
        assert response.status_code == 400
        assert "blocked" in response.json()["detail"].lower()

    def test_scrape_rejects_file_scheme(self, client):
        """file:// scheme should be blocked."""
        response = client.post(
            "/api/extraction/scrape",
            json={"url": "file:///etc/passwd"},
        )
        assert response.status_code == 400
        assert "blocked" in response.json()["detail"].lower()

    def test_scrape_rejects_ipv6_loopback(self, client):
        """IPv6 loopback should be blocked."""
        response = client.post(
            "/api/extraction/scrape",
            json={"url": "http://[::1]/"},
        )
        assert response.status_code == 400
        assert "blocked" in response.json()["detail"].lower()


class TestOCRLanguageValidation:
    """OCR language validation tests."""

    def test_ocr_rejects_invalid_language(self, client):
        """Invalid language code should be rejected."""
        response = client.post(
            "/api/extraction/ocr?language=xyz",
            files={"file": ("test.png", b"fake", "image/png")},
        )
        assert response.status_code == 400
        assert "unsupported language" in response.json()["detail"].lower()