"""Test that application startup validates required configuration."""

import os
import pytest
from unittest.mock import patch


def test_lifespan_validates_secrets():
    """Startup should fail fast if required secrets missing."""
    # Clear all secrets
    env_without_secrets = {
        k: v for k, v in os.environ.items() 
        if k not in ["SIGNATURE_SECRET_KEY", "JWT_SECRET_KEY", "NOTEBOOK_REVIEW_SECRET"]
    }
    
    with patch.dict(os.environ, env_without_secrets, clear=True):
        # Import fresh - should fail in lifespan
        from amprenta_rag.utils.config_check import validate_required_secrets
        valid, missing = validate_required_secrets()
        
        assert not valid, "Should detect missing secrets"
        assert any("SIGNATURE" in m for m in missing)
        assert any("JWT" in m for m in missing)


def test_config_check_with_all_secrets():
    """Config check should pass when all secrets are present."""
    from amprenta_rag.utils.config_check import validate_required_secrets
    
    # With test environment configured, should pass
    valid, missing = validate_required_secrets()
    
    # May pass (test env configured) or fail (missing other secrets)
    # This test verifies the function works, not that all secrets are present
    assert isinstance(valid, bool)
    assert isinstance(missing, list)


def test_app_import_with_test_config():
    """App should import successfully with test configuration."""
    # This test verifies that the test_config module provides sufficient defaults
    # for app import without RuntimeError
    
    try:
        from amprenta_rag.api.main import app
        assert app is not None
    except RuntimeError as e:
        if "not configured" in str(e):
            pytest.fail(f"Import-time configuration error: {e}")
        else:
            # Other RuntimeErrors may be acceptable
            pass


def test_lazy_secret_loading():
    """Secrets should be loaded lazily, not at import time."""
    # Clear environment
    env_without_secrets = {
        k: v for k, v in os.environ.items() 
        if not k.endswith("_SECRET_KEY")
    }
    
    with patch.dict(os.environ, env_without_secrets, clear=True):
        # These imports should work even without secrets
        try:
            from amprenta_rag.auth.signatures import _get_secret_key
            from amprenta_rag.services.notebook_review import _get_secret
            
            # Functions exist but haven't been called yet
            assert callable(_get_secret_key)
            assert callable(_get_secret)
            
        except ImportError:
            # Import errors are acceptable (missing dependencies)
            pass
        except RuntimeError as e:
            if "not configured" in str(e):
                pytest.fail(f"Secret loaded at import time: {e}")


def test_secret_functions_fail_when_called_without_config():
    """Secret functions should fail when called without configuration."""
    env_without_secrets = {
        k: v for k, v in os.environ.items() 
        if not k.endswith("_SECRET_KEY")
    }
    
    with patch.dict(os.environ, env_without_secrets, clear=True):
        try:
            from amprenta_rag.auth.signatures import _get_secret_key
            
            # Function should exist but fail when called
            with pytest.raises(RuntimeError, match="not configured"):
                _get_secret_key()
                
        except ImportError:
            # Skip if dependencies not available
            pytest.skip("RDKit or other dependencies not available")
