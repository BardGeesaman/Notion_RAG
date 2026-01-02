"""Tests for input validation and sanitization."""
import pytest


class TestHTMLSanitization:
    """Test HTML sanitization utilities."""
    
    def test_escape_html_removes_script_tags(self):
        """Test that script tags are escaped."""
        from amprenta_rag.utils.sanitize import escape_html
        
        malicious = '<script>alert("xss")</script>Hello'
        result = escape_html(malicious)
        
        assert '<script>' not in result
        assert 'alert' not in result or '&lt;script&gt;' in result
        assert 'Hello' in result
    
    def test_escape_html_removes_onclick(self):
        """Test that event handlers are removed."""
        from amprenta_rag.utils.sanitize import escape_html
        
        malicious = '<div onclick="evil()">Click me</div>'
        result = escape_html(malicious)
        
        assert 'onclick' not in result
        assert 'Click me' in result
    
    def test_sanitize_html_allows_safe_tags(self):
        """Test that safe tags are preserved."""
        from amprenta_rag.utils.sanitize import sanitize_html
        
        safe = '<b>Bold</b> and <i>italic</i>'
        result = sanitize_html(safe)
        
        assert '<b>' in result
        assert '<i>' in result
    
    def test_sanitize_html_removes_unsafe_tags(self):
        """Test that unsafe tags are removed."""
        from amprenta_rag.utils.sanitize import sanitize_html
        
        unsafe = '<script>bad</script><b>good</b>'
        result = sanitize_html(unsafe)
        
        assert '<script>' not in result
        assert '<b>good</b>' in result
    
    def test_highlight_mentions_prevents_xss(self):
        """Test that _highlight_mentions prevents XSS."""
        from scripts.dashboard.components.comment_widget import _highlight_mentions
        
        malicious = '<script>alert("xss")</script>@user hello'
        result = _highlight_mentions(malicious)
        
        assert '<script>' not in result
        # @user should still be highlighted
        assert '@user' in result or 'user' in result


class TestPydanticStrictMode:
    """Test Pydantic strict mode validation."""
    
    def test_strict_schema_rejects_type_coercion(self):
        """Test that strict schemas don't coerce types."""
        from amprenta_rag.api.schemas import StrictBaseSchema
        from pydantic import Field
        
        class TestStrict(StrictBaseSchema):
            count: int
        
        # Strict mode should reject string "1" for int field
        with pytest.raises(Exception):  # ValidationError
            TestStrict(count="1")
    
    def test_validate_smiles_rejects_invalid(self):
        """Test SMILES validation rejects invalid input."""
        from amprenta_rag.api.schemas import validate_smiles
        
        with pytest.raises(ValueError):
            validate_smiles("")
        
        with pytest.raises(ValueError):
            validate_smiles("<script>alert('xss')</script>")
    
    def test_validate_url_blocks_javascript(self):
        """Test URL validation blocks dangerous schemes."""
        from amprenta_rag.api.schemas import validate_url
        
        with pytest.raises(ValueError):
            validate_url("javascript:alert('xss')")
        
        # Safe URLs should pass
        assert validate_url("https://example.com") == "https://example.com"


class TestRequestSizeLimits:
    """Test request body size limits."""
    
    def test_size_limit_middleware_configured(self):
        """Test that size limit middleware is configured."""
        from amprenta_rag.api.main import app
        
        middleware_classes = [m.cls.__name__ for m in app.user_middleware]
        assert "RequestSizeLimitMiddleware" in middleware_classes
    
    def test_upload_routes_have_higher_limit(self):
        """Test that upload routes are configured for higher limits."""
        from amprenta_rag.api.middleware import UPLOAD_ROUTES, MAX_UPLOAD_SIZE, MAX_REQUEST_SIZE
        
        # Upload limit should be higher than default
        assert MAX_UPLOAD_SIZE > MAX_REQUEST_SIZE
        
        # Should have upload routes defined
        assert len(UPLOAD_ROUTES) > 0
        assert "/api/v1/datasets/upload" in UPLOAD_ROUTES or any(
            "upload" in route for route in UPLOAD_ROUTES
        )
    
    def test_size_limit_returns_413(self):
        """Test that oversized requests return 413."""
        from amprenta_rag.api.middleware import RequestSizeLimitMiddleware
        from starlette.testclient import TestClient
        from starlette.applications import Starlette
        from starlette.responses import PlainTextResponse
        from starlette.routing import Route
        
        async def echo(request):
            body = await request.body()
            return PlainTextResponse(f"Received {len(body)} bytes")
        
        test_app = Starlette(routes=[Route("/test", echo, methods=["POST"])])
        test_app.add_middleware(RequestSizeLimitMiddleware, default_limit=100)  # 100 bytes
        
        client = TestClient(test_app)
        
        # Small request should work
        response = client.post("/test", content="small")
        assert response.status_code == 200
        
        # Large request should be rejected
        response = client.post("/test", content="x" * 200, headers={"content-length": "200"})
        assert response.status_code == 413
        assert "request_too_large" in response.json()["error"]
