"""Tests for electronic signature service."""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest

from amprenta_rag.auth.signatures import (
    create_signature,
    get_signatures,
    verify_signature,
)


class TestCreateSignature:
    """Tests for create_signature."""

    @patch("amprenta_rag.auth.signatures.verify_password")
    def test_create_signature_success(self, mock_verify):
        """Test successful signature creation."""
        mock_verify.return_value = True
        
        user_id = uuid4()
        entity_id = uuid4()
        
        mock_user = MagicMock()
        mock_user.id = user_id
        mock_user.username = "testuser"
        mock_user.password_hash = "hashed_password"
        
        mock_signature = MagicMock()
        mock_signature.id = uuid4()
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_user
        
        result = create_signature(
            str(user_id),
            "approve",
            "compound",
            str(entity_id),
            "correct_password",
            "I approve this compound",
            mock_db,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called

    @patch("amprenta_rag.auth.signatures.verify_password")
    def test_create_signature_wrong_password(self, mock_verify):
        """Test signature creation with wrong password."""
        mock_verify.return_value = False
        
        user_id = uuid4()
        
        mock_user = MagicMock()
        mock_user.password_hash = "hashed_password"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_user
        
        result = create_signature(
            str(user_id),
            "approve",
            "compound",
            str(uuid4()),
            "wrong_password",
            "I approve",
            mock_db,
        )
        
        assert result is None


class TestVerifySignature:
    """Tests for verify_signature."""

    @patch("amprenta_rag.auth.signatures.hmac.new")
    def test_verify_signature_valid(self, mock_hmac):
        """Test verification of valid signature."""
        mock_hash = MagicMock()
        mock_hash.hexdigest.return_value = "abc123def456"
        mock_hmac.return_value = mock_hash
        
        mock_signature = MagicMock()
        mock_signature.id = uuid4()
        mock_signature.user_id = uuid4()
        mock_signature.signature_hash = "abc123def456"
        mock_signature.action = "approve"
        mock_signature.entity_type = "compound"
        mock_signature.entity_id = uuid4()
        mock_signature.meaning = "I approve"
        mock_signature.verified_at = None
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_signature
        
        result = verify_signature(str(mock_signature.id), mock_db)
        
        assert result is True

    @patch("amprenta_rag.auth.signatures.hmac.new")
    def test_verify_signature_tampered(self, mock_hmac):
        """Test verification of tampered signature."""
        mock_hash = MagicMock()
        mock_hash.hexdigest.return_value = "expected_hash"
        mock_hmac.return_value = mock_hash
        
        mock_signature = MagicMock()
        mock_signature.signature_hash = "tampered_hash"
        mock_signature.user_id = uuid4()
        mock_signature.action = "approve"
        mock_signature.entity_type = "compound"
        mock_signature.entity_id = uuid4()
        mock_signature.meaning = "I approve"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_signature
        
        result = verify_signature(str(uuid4()), mock_db)
        
        assert result is False


class TestGetSignatures:
    """Tests for get_signatures."""

    def test_get_signatures_returns_list(self):
        """Test that get_signatures returns list."""
        mock_sig = MagicMock()
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = [mock_sig]
        
        result = get_signatures("compound", str(uuid4()), mock_db)
        
        assert isinstance(result, list)
        assert len(result) == 1

    def test_get_signatures_empty(self):
        """Test get_signatures with no signatures."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = []
        
        result = get_signatures("compound", str(uuid4()), mock_db)
        
        assert result == []

