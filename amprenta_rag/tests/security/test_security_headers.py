"""Tests for security headers middleware."""
from fastapi.testclient import TestClient
from amprenta_rag.api.main import app


def test_x_content_type_options_header():
    """Test X-Content-Type-Options header is present."""
    client = TestClient(app)
    response = client.get("/health")
    assert response.headers.get("X-Content-Type-Options") == "nosniff"


def test_x_frame_options_header():
    """Test X-Frame-Options header is present."""
    client = TestClient(app)
    response = client.get("/health")
    assert response.headers.get("X-Frame-Options") == "SAMEORIGIN"


def test_referrer_policy_header():
    """Test Referrer-Policy header is present."""
    client = TestClient(app)
    response = client.get("/health")
    assert response.headers.get("Referrer-Policy") == "strict-origin-when-cross-origin"


def test_csp_header():
    """Test Content-Security-Policy header is present."""
    client = TestClient(app)
    response = client.get("/health")
    csp = response.headers.get("Content-Security-Policy")
    assert csp is not None
    assert "default-src 'self'" in csp
