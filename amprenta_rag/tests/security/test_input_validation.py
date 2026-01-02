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
