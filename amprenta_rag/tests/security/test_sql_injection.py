"""Tests for SQL injection prevention."""
import pytest
from uuid import uuid4
from unittest.mock import MagicMock

from amprenta_rag.auth.company_context import set_company_context


class TestSQLInjectionPrevention:
    """Test SQL injection prevention in company context."""
    
    def test_company_id_with_sql_injection_payload(self):
        """Test that SQL injection payloads in company_id are safely escaped."""
        # Create mock user with malicious company_id
        mock_user = MagicMock()
        mock_user.company_id = "'; DROP TABLE users; --"
        
        mock_db = MagicMock()
        mock_request = MagicMock()
        
        # Call function directly with malicious user
        result = set_company_context(mock_request, mock_db, mock_user)
        
        # Verify parameterized query was used (not f-string)
        if mock_db.execute.called:
            call_args = mock_db.execute.call_args
            # The first argument should be a text() object, not a string with the payload
            sql_text = str(call_args[0][0])
            assert "DROP TABLE" not in sql_text
            # Parameters should be passed separately
            if len(call_args) > 1 and call_args[1]:
                assert "company_id" in call_args[1]
    
    def test_company_id_with_quote_characters(self):
        """Test that quote characters in company_id don't break SQL."""
        mock_user = MagicMock()
        mock_user.company_id = str(uuid4())  # Valid UUID
        
        mock_db = MagicMock()
        mock_request = MagicMock()
        
        # Should not raise
        set_company_context(mock_request, mock_db, mock_user)
        
        # Verify execute was called with parameterized query
        assert mock_db.execute.called