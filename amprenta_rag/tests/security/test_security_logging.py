"""A09: Security Logging tests."""

import pytest
import logging
from unittest.mock import patch, MagicMock
from uuid import uuid4


class TestSecurityLogger:
    """Test security logging functionality."""
    
    def test_security_logger_exists(self):
        """Security logger module should exist."""
        from amprenta_rag.auth.security_logger import security_logger
        
        assert security_logger is not None
        assert security_logger.name == "security"
    
    def test_log_auth_failure(self):
        """Auth failure should be logged."""
        from amprenta_rag.auth.security_logger import log_auth_failure
        
        with patch('amprenta_rag.auth.security_logger.security_logger') as mock_logger:
            log_auth_failure(
                username="test_user",
                ip_address="192.168.1.1",
                reason="invalid_password"
            )
            
            mock_logger.log.assert_called_once()
            call_args = mock_logger.log.call_args
            
            # Check log level is WARNING
            assert call_args[0][0] == logging.WARNING
            
            # Check message contains key info
            message = call_args[0][1]
            assert "auth_failure" in message
            assert "192.168.1.1" in message
    
    def test_log_auth_success(self):
        """Auth success should be logged."""
        from amprenta_rag.auth.security_logger import log_auth_success
        
        user_id = uuid4()
        
        with patch('amprenta_rag.auth.security_logger.security_logger') as mock_logger:
            log_auth_success(user_id=user_id, ip_address="10.0.0.1")
            
            mock_logger.log.assert_called_once()
            message = mock_logger.log.call_args[0][1]
            assert "auth_success" in message
            assert str(user_id) in message
    
    def test_log_rate_limit(self):
        """Rate limit events should be logged."""
        from amprenta_rag.auth.security_logger import log_rate_limit
        
        with patch('amprenta_rag.auth.security_logger.security_logger') as mock_logger:
            log_rate_limit(
                user_id=uuid4(),
                endpoint="/api/v1/signatures",
                ip_address="172.16.0.1",
                limit="5/minute"
            )
            
            mock_logger.log.assert_called_once()
            message = mock_logger.log.call_args[0][1]
            assert "rate_limit" in message
    
    def test_log_account_locked(self):
        """Account lockout should be logged."""
        from amprenta_rag.auth.security_logger import log_account_locked
        
        user_id = uuid4()
        
        with patch('amprenta_rag.auth.security_logger.security_logger') as mock_logger:
            log_account_locked(
                user_id=user_id,
                ip_address="192.168.1.100",
                failed_attempts=5
            )
            
            mock_logger.log.assert_called_once()
            message = mock_logger.log.call_args[0][1]
            assert "account_locked" in message
            assert str(user_id) in message
    
    def test_security_event_enum(self):
        """Security event types should be defined."""
        from amprenta_rag.auth.security_logger import SecurityEvent
        
        # Verify key events exist
        assert SecurityEvent.AUTH_SUCCESS.value == "auth_success"
        assert SecurityEvent.AUTH_FAILURE.value == "auth_failure"
        assert SecurityEvent.AUTHZ_FAILURE.value == "authz_failure"
        assert SecurityEvent.RATE_LIMIT_HIT.value == "rate_limit_hit"
        assert SecurityEvent.ACCOUNT_LOCKED.value == "account_locked"
