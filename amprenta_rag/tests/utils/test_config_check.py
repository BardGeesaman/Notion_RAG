"""Tests for configuration validation utilities."""

import os
from unittest.mock import patch
from amprenta_rag.utils.config_check import validate_required_secrets


class TestValidateRequiredSecrets:
    """Tests for secret validation function."""
    
    def test_all_secrets_present(self):
        """Should return True when all required secrets present."""
        with patch.dict(os.environ, {
            "SIGNATURE_SECRET_KEY": "test-key",
            "JWT_SECRET_KEY": "test-jwt",
            "POSTGRES_PASSWORD": "test-pass",
        }):
            valid, missing = validate_required_secrets()
            assert valid is True
            assert missing == []
    
    def test_missing_signature_key(self):
        """Should return False when SIGNATURE_SECRET_KEY missing."""
        with patch.dict(os.environ, {
            "JWT_SECRET_KEY": "test-jwt",
            "POSTGRES_PASSWORD": "test-pass",
        }, clear=True):
            valid, missing = validate_required_secrets()
            assert valid is False
            assert any("SIGNATURE_SECRET_KEY" in item for item in missing)
    
    def test_missing_jwt_key(self):
        """Should return False when JWT_SECRET_KEY missing."""
        with patch.dict(os.environ, {
            "SIGNATURE_SECRET_KEY": "test-key",
            "POSTGRES_PASSWORD": "test-pass",
        }, clear=True):
            valid, missing = validate_required_secrets()
            assert valid is False
            assert any("JWT_SECRET_KEY" in item for item in missing)
    
    def test_postgres_url_overrides_password(self):
        """Should not require POSTGRES_PASSWORD when POSTGRES_URL present."""
        with patch.dict(os.environ, {
            "SIGNATURE_SECRET_KEY": "test-key",
            "JWT_SECRET_KEY": "test-jwt",
            "POSTGRES_URL": "postgresql://user:pass@host/db",
        }, clear=True):
            valid, missing = validate_required_secrets()
            assert valid is True
            assert missing == []
    
    def test_database_url_overrides_password(self):
        """Should not require POSTGRES_PASSWORD when DATABASE_URL present."""
        with patch.dict(os.environ, {
            "SIGNATURE_SECRET_KEY": "test-key",
            "JWT_SECRET_KEY": "test-jwt",
            "DATABASE_URL": "postgresql://user:pass@host/db",
        }, clear=True):
            valid, missing = validate_required_secrets()
            assert valid is True
            assert missing == []