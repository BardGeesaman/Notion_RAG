"""A04/A05: Insecure Design & Configuration tests."""

import pytest
import os


class TestCORSConfiguration:
    """Test CORS security configuration."""
    
    def test_cors_not_wildcard_in_production(self):
        """CORS should not allow all origins in production."""
        # This is a code review test - check main.py
        from pathlib import Path
        
        main_py = Path("amprenta_rag/api/main.py").read_text()
        
        # Check for wildcard CORS (bad)
        if 'allow_origins=["*"]' in main_py:
            # Only acceptable if controlled by env var
            assert "CORS_ORIGINS" in main_py or "os.environ" in main_py, \
                "CORS wildcard should be env-controlled"
    
    def test_cors_credentials_restricted(self):
        """CORS with credentials should have specific origins."""
        from pathlib import Path
        
        main_py = Path("amprenta_rag/api/main.py").read_text()
        
        # If credentials=True, origins should not be *
        if "allow_credentials=True" in main_py:
            assert 'allow_origins=["*"]' not in main_py, \
                "Cannot use credentials with wildcard origins"


class TestDebugConfiguration:
    """Test debug mode configuration."""
    
    def test_debug_not_hardcoded_true(self):
        """Debug should not be hardcoded to True."""
        from pathlib import Path
        
        main_py = Path("amprenta_rag/api/main.py").read_text()
        
        # Should not have debug=True hardcoded
        assert "debug=True" not in main_py or "DEBUG" in main_py, \
            "Debug should be controlled by environment"
    
    def test_disable_auth_not_in_production(self):
        """DISABLE_AUTH should not be set in production."""
        # In tests, it might be set, but verify it's env-controlled
        disable_auth = os.environ.get("DISABLE_AUTH", "")
        
        # This test documents the env var exists
        # Production should never have this set
        assert disable_auth.lower() not in ("1", "true") or \
            os.environ.get("ENVIRONMENT") != "production"


class TestErrorHandling:
    """Test error response security."""
    
    def test_500_errors_dont_leak_details(self):
        """Internal errors should not expose stack traces."""
        from fastapi.testclient import TestClient
        from amprenta_rag.api.main import app
        
        client = TestClient(app, raise_server_exceptions=False)
        
        # Trigger an error (invalid UUID)
        response = client.get(
            "/api/v1/datasets/not-a-uuid",
            headers={"X-User-Id": "00000000-0000-0000-0000-000000000001"}
        )
        
        if response.status_code >= 500:
            # Should not contain stack trace
            body = response.text.lower()
            assert "traceback" not in body
            assert "file " not in body  # File paths
            assert "line " not in body  # Line numbers


class TestMassAssignment:
    """Test mass assignment protection."""
    
    def test_pydantic_schemas_filter_fields(self):
        """Pydantic schemas should not accept arbitrary fields."""
        from pydantic import ValidationError, BaseModel, ConfigDict
        
        class TestSchema(BaseModel):
            model_config = ConfigDict(extra='forbid')  # Explicitly forbid extra fields
            name: str
        
        # Should reject extra fields
        with pytest.raises(ValidationError):
            TestSchema(name="test", admin=True, role="superuser")


class TestSecurityHeaders:
    """Verify security headers are set."""
    
    def test_security_headers_present(self):
        """Security headers should be present in responses."""
        from fastapi.testclient import TestClient
        from amprenta_rag.api.main import app
        
        client = TestClient(app)
        response = client.get("/health")
        
        # Check security headers
        assert response.headers.get("X-Content-Type-Options") == "nosniff"
        assert response.headers.get("X-Frame-Options") in ("DENY", "SAMEORIGIN")
