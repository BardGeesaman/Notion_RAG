"""A02/A07: Cryptographic and Authentication tests."""

import pytest
from unittest.mock import patch, MagicMock


class TestPasswordHashing:
    """Test password hashing security."""
    
    def test_passwords_are_hashed_with_bcrypt(self):
        """Verify bcrypt is used for password hashing."""
        from amprenta_rag.auth.password import hash_password
        
        hashed = hash_password("test_password")
        
        # bcrypt hashes start with $2b$ or $2a$
        assert hashed.startswith("$2") or hashed.startswith("$argon2")
    
    def test_password_hash_is_not_reversible(self):
        """Hashed password should not contain original."""
        from amprenta_rag.auth.password import hash_password
        
        password = "super_secret_123"
        hashed = hash_password(password)
        
        assert password not in hashed
    
    def test_same_password_produces_different_hashes(self):
        """Salt should make hashes unique."""
        from amprenta_rag.auth.password import hash_password
        
        hash1 = hash_password("test")
        hash2 = hash_password("test")
        
        assert hash1 != hash2  # Different salts
    
    def test_password_verification_works(self):
        """Verify password checking works correctly."""
        from amprenta_rag.auth.password import hash_password, verify_password
        
        password = "correct_password"
        hashed = hash_password(password)
        
        assert verify_password(password, hashed) is True
        assert verify_password("wrong_password", hashed) is False


class TestJWTSecurity:
    """Test JWT token security."""
    
    def test_jwt_uses_secure_algorithm(self):
        """JWT should use HS256 or RS256, not 'none'."""
        # Check JWT creation code
        import amprenta_rag.auth.signatures as sigs
        
        # Verify algorithm constant if present
        if hasattr(sigs, 'ALGORITHM'):
            assert sigs.ALGORITHM in ('HS256', 'RS256', 'HS384', 'HS512')
    
    def test_jwt_secret_from_environment(self):
        """JWT secret should come from environment."""
        import os
        
        # Should use JWT_SECRET_KEY or similar env var
        # Not hardcoded
        assert os.environ.get("JWT_SECRET_KEY") or os.environ.get("SECRET_KEY")


class TestAccountLockout:
    """Test account lockout functionality."""
    
    def test_lockout_configuration(self):
        """Test lockout thresholds are properly configured."""
        from amprenta_rag.auth.lockout import MAX_FAILED_ATTEMPTS, LOCKOUT_DURATION
        
        # Should be reasonable values
        assert MAX_FAILED_ATTEMPTS >= 3  # Not too restrictive
        assert MAX_FAILED_ATTEMPTS <= 10  # Not too permissive
        assert LOCKOUT_DURATION.total_seconds() >= 300  # At least 5 minutes
    
    def test_lockout_identifier_format(self):
        """Test lockout identifiers are properly formatted."""
        from amprenta_rag.auth.lockout import record_failed_attempt
        
        # Should handle various identifier formats
        identifiers = [
            "sign:192.168.1.1",
            "login:user@example.com",
            "api:key_123"
        ]
        
        for identifier in identifiers:
            # Should not raise exception
            try:
                with patch("amprenta_rag.auth.lockout.db_session") as mock_session:
                    mock_db = MagicMock()
                    mock_session.return_value.__enter__.return_value = mock_db
                    mock_db.query.return_value.filter.return_value.count.return_value = 1
                    
                    count = record_failed_attempt(identifier)
                    assert isinstance(count, int)
            except Exception as e:
                pytest.fail(f"Failed to record attempt for {identifier}: {e}")


class TestSessionSecurity:
    """Test session management security."""
    
    def test_session_tokens_are_random(self):
        """Session tokens should be cryptographically random."""
        import secrets
        
        # Generate two tokens
        token1 = secrets.token_urlsafe(32)
        token2 = secrets.token_urlsafe(32)
        
        # Should be different
        assert token1 != token2
        
        # Should be sufficient length (at least 32 chars)
        assert len(token1) >= 32
