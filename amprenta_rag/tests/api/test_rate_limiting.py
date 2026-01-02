"""Tests for rate limiting infrastructure."""
import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.api.main import app


def test_rate_limiter_configured():
    """Test that rate limiter is configured on app."""
    assert hasattr(app.state, "limiter")


def test_rate_limit_exceeded_handler():
    """Test rate limit exceeded returns 429."""
    from amprenta_rag.api.rate_limit import rate_limit_exceeded_handler
    from fastapi import Request
    
    request = MagicMock(spec=Request)
    
    # Create mock exception with required attributes
    exc = MagicMock()
    exc.detail = "5 per minute"
    exc.retry_after = 60
    
    response = rate_limit_exceeded_handler(request, exc)
    
    assert response.status_code == 429
    assert "Retry-After" in response.headers


def test_get_user_or_ip_with_user():
    """Test get_user_or_ip returns user ID when authenticated."""
    from amprenta_rag.api.rate_limit import get_user_or_ip, set_user_state
    
    request = MagicMock()
    request.state = MagicMock()
    
    mock_user = MagicMock()
    mock_user.id = uuid4()
    set_user_state(request, mock_user)
    
    result = get_user_or_ip(request)
    assert result == f"user:{mock_user.id}"


def test_get_user_or_ip_without_user():
    """Test get_user_or_ip returns IP when not authenticated."""
    from amprenta_rag.api.rate_limit import get_user_or_ip
    
    request = MagicMock()
    request.state = MagicMock(spec=[])  # No user attribute
    request.client.host = "192.168.1.1"
    
    with patch("amprenta_rag.api.rate_limit.get_remote_address", return_value="192.168.1.1"):
        result = get_user_or_ip(request)
    
    assert result == "ip:192.168.1.1"


def test_signature_endpoint_rate_limited():
    """Test that signature endpoint has rate limiting."""
    from amprenta_rag.api.routers.signatures import sign_document
    
    # Check the endpoint has the limiter decorator
    assert hasattr(sign_document, "__wrapped__") or hasattr(sign_document, "_limit")


def test_share_link_validate_rate_limited():
    """Test that share link validation has rate limiting."""
    from amprenta_rag.api.routers.share_links import validate_token
    
    assert hasattr(validate_token, "__wrapped__") or hasattr(validate_token, "_limit")


def test_timeout_middleware_configured():
    """Test that timeout middleware is configured on app."""
    from amprenta_rag.api.main import app
    
    # Check middleware stack includes TimeoutMiddleware
    middleware_classes = [m.cls.__name__ for m in app.user_middleware]
    assert "TimeoutMiddleware" in middleware_classes


def test_timeout_middleware_returns_504_on_timeout():
    """Test timeout middleware returns 504 on timeout."""
    import asyncio
    from amprenta_rag.api.middleware import TimeoutMiddleware
    from starlette.testclient import TestClient
    from starlette.applications import Starlette
    from starlette.responses import PlainTextResponse
    from starlette.routing import Route
    
    async def slow_endpoint(request):
        await asyncio.sleep(5)  # Simulate slow response
        return PlainTextResponse("done")
    
    test_app = Starlette(routes=[Route("/slow", slow_endpoint)])
    test_app.add_middleware(TimeoutMiddleware, timeout=1)  # 1 second timeout
    
    client = TestClient(test_app)
    response = client.get("/slow")
    
    assert response.status_code == 504
    assert "request_timeout" in response.json()["error"]
