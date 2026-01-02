"""Tests for secrets management utilities."""

import os
from unittest.mock import patch

import pytest

from amprenta_rag.utils.secrets import (
    get_auth_secret,
    get_email_credential,
    clear_cache,
    get_environment_info,
)


class TestGetAuthSecret:
    """Tests for get_auth_secret function."""

    def test_get_signature_key_from_env(self):
        """Test retrieving signature key from environment variable."""
        with patch.dict(os.environ, {"SIGNATURE_SECRET_KEY": "test-sig-key-123"}):
            result = get_auth_secret("signature_key")
            assert result == "test-sig-key-123"

    def test_get_jwt_key_from_env(self):
        """Test retrieving JWT key from environment variable."""
        with patch.dict(os.environ, {"JWT_SECRET_KEY": "test-jwt-key-456"}):
            result = get_auth_secret("jwt_key")
            assert result == "test-jwt-key-456"

    def test_unknown_secret_returns_none(self):
        """Test that unknown secret names return None."""
        result = get_auth_secret("unknown_secret")
        assert result is None

    def test_missing_env_var_returns_none(self):
        """Test that missing environment variable returns None."""
        with patch.dict(os.environ, {}, clear=True):
            # Ensure the env var is not set
            os.environ.pop("SIGNATURE_SECRET_KEY", None)
            result = get_auth_secret("signature_key")
            assert result is None


class TestGetEmailCredential:
    """Tests for get_email_credential function."""

    def test_get_smtp_user_from_env(self):
        """Test retrieving SMTP user from environment variable."""
        with patch.dict(os.environ, {"SMTP_USER": "test@example.com"}):
            result = get_email_credential("smtp_user")
            assert result == "test@example.com"

    def test_get_smtp_password_from_env(self):
        """Test retrieving SMTP password from environment variable."""
        with patch.dict(os.environ, {"SMTP_PASSWORD": "secret-password"}):
            result = get_email_credential("smtp_password")
            assert result == "secret-password"

    def test_get_smtp_host_from_env(self):
        """Test retrieving SMTP host from environment variable."""
        with patch.dict(os.environ, {"SMTP_HOST": "mail.example.com"}):
            result = get_email_credential("smtp_host")
            assert result == "mail.example.com"

    def test_get_smtp_port_from_env(self):
        """Test retrieving SMTP port from environment variable."""
        with patch.dict(os.environ, {"SMTP_PORT": "465"}):
            result = get_email_credential("smtp_port")
            assert result == "465"

    def test_get_from_email_from_env(self):
        """Test retrieving from email from environment variable."""
        with patch.dict(os.environ, {"FROM_EMAIL": "noreply@example.com"}):
            result = get_email_credential("from_email")
            assert result == "noreply@example.com"

    def test_unknown_credential_returns_none(self):
        """Test that unknown credential names return None."""
        result = get_email_credential("unknown_credential")
        assert result is None


class TestSecurityFailFast:
    """Tests to verify fail-fast behavior for missing secrets."""

    def test_signatures_module_requires_secret(self):
        """Test that signatures module fails without SIGNATURE_SECRET_KEY."""
        # Clear environment and reload module
        with patch.dict(os.environ, {}, clear=True):
            os.environ.pop("SIGNATURE_SECRET_KEY", None)
            
            # Importing should raise RuntimeError
            import importlib
            import sys
            
            # Remove cached module
            if "amprenta_rag.auth.signatures" in sys.modules:
                del sys.modules["amprenta_rag.auth.signatures"]
            
            with pytest.raises(RuntimeError, match="SIGNATURE_SECRET_KEY not configured"):
                importlib.import_module("amprenta_rag.auth.signatures")


class TestClearCache:
    """Tests for cache management."""

    def test_clear_cache_resets_state(self):
        """Test that clear_cache resets the secrets cache."""
        # This should not raise
        clear_cache()
        
        # Verify cache is empty via environment info
        info = get_environment_info()
        assert info["cache_size"] == 0


class TestGetEnvironmentInfo:
    """Tests for environment info function."""

    def test_returns_expected_keys(self):
        """Test that get_environment_info returns expected keys."""
        info = get_environment_info()
        
        assert "is_aws" in info
        assert "environment" in info
        assert "cache_size" in info
        assert isinstance(info["is_aws"], bool)
        assert isinstance(info["cache_size"], int)
