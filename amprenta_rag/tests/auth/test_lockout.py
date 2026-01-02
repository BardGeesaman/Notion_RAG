"""Tests for account lockout functionality."""
import pytest
from datetime import datetime, timedelta, timezone
from unittest.mock import patch, MagicMock
from uuid import uuid4

from amprenta_rag.auth.lockout import (
    record_failed_attempt,
    is_locked_out,
    clear_failed_attempts,
    MAX_FAILED_ATTEMPTS,
    LOCKOUT_DURATION,
)


class TestLockout:
    """Test lockout functionality."""
    
    def test_record_failed_attempt(self):
        """Test recording a failed attempt."""
        identifier = f"test:{uuid4().hex[:8]}"
        
        with patch("amprenta_rag.auth.lockout.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.count.return_value = 1
            
            count = record_failed_attempt(identifier)
            
            assert mock_db.add.called
            assert mock_db.commit.called
    
    def test_is_locked_out_false_when_no_attempts(self):
        """Test not locked when no failed attempts."""
        identifier = f"test:{uuid4().hex[:8]}"
        
        with patch("amprenta_rag.auth.lockout.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.count.return_value = 0
            
            locked, remaining = is_locked_out(identifier)
            
            assert locked is False
            assert remaining is None
    
    def test_is_locked_out_true_after_max_attempts(self):
        """Test locked after max failed attempts."""
        identifier = f"test:{uuid4().hex[:8]}"
        
        with patch("amprenta_rag.auth.lockout.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.count.return_value = MAX_FAILED_ATTEMPTS
            
            # Mock last attempt time as recent
            recent_time = datetime.now(timezone.utc) - timedelta(minutes=1)
            mock_db.query.return_value.filter.return_value.scalar.return_value = recent_time
            
            locked, remaining = is_locked_out(identifier)
            
            assert locked is True
            assert remaining is not None
            assert remaining > 0
    
    def test_clear_failed_attempts(self):
        """Test clearing failed attempts."""
        identifier = f"test:{uuid4().hex[:8]}"
        
        with patch("amprenta_rag.auth.lockout.db_session") as mock_session:
            mock_db = MagicMock()
            mock_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.delete.return_value = 3
            
            clear_failed_attempts(identifier)
            
            assert mock_db.query.return_value.filter.return_value.delete.called
            assert mock_db.commit.called
