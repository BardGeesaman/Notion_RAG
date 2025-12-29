"""Tests for electronic signatures API endpoints."""

from __future__ import annotations

from datetime import datetime
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestSignDocument:
    """Tests for POST /api/v1/signatures/sign endpoint."""

    @patch("amprenta_rag.api.routers.signatures.create_signature")
    def test_sign_success(self, mock_create):
        """Test successful signature creation."""
        mock_sig = MagicMock()
        mock_sig.id = uuid4()
        mock_sig.signature_hash = "abc123def456"
        mock_sig.timestamp = datetime(2024, 1, 1)
        
        mock_create.return_value = mock_sig
        
        response = client.post(
            "/api/v1/signatures/sign",
            json={
                "user_id": str(uuid4()),
                "action": "approve",
                "entity_type": "compound",
                "entity_id": str(uuid4()),
                "meaning": "I approve this compound",
                "password": "correct_password",
            },
        )
        
        assert response.status_code == 201
        data = response.json()
        assert data["valid"] is True
        assert "signature_hash" in data

    @patch("amprenta_rag.api.routers.signatures.create_signature")
    def test_sign_wrong_password(self, mock_create):
        """Test signature creation with wrong password."""
        mock_create.return_value = None  # Password verification failed
        
        response = client.post(
            "/api/v1/signatures/sign",
            json={
                "user_id": str(uuid4()),
                "action": "approve",
                "entity_type": "compound",
                "entity_id": str(uuid4()),
                "meaning": "I approve",
                "password": "wrong_password",
            },
        )
        
        assert response.status_code == 401
        assert "authentication failed" in response.json()["detail"].lower()


class TestVerifySignature:
    """Tests for GET /api/v1/signatures/{id}/verify endpoint."""

    @patch("amprenta_rag.api.routers.signatures.verify_signature")
    def test_verify_valid_signature(self, mock_verify):
        """Test verification of valid signature."""
        mock_verify.return_value = True
        
        signature_id = uuid4()
        response = client.get(f"/api/v1/signatures/{signature_id}/verify")
        
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True

    @patch("amprenta_rag.api.routers.signatures.verify_signature")
    def test_verify_invalid_signature(self, mock_verify):
        """Test verification of invalid signature."""
        mock_verify.return_value = False
        
        signature_id = uuid4()
        response = client.get(f"/api/v1/signatures/{signature_id}/verify")
        
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False


class TestListSignatures:
    """Tests for GET /api/v1/signatures/{entity_type}/{entity_id} endpoint."""

    @patch("amprenta_rag.api.routers.signatures.get_signatures")
    def test_list_signatures_for_entity(self, mock_get_sigs):
        """Test listing signatures for entity."""
        mock_sig = MagicMock()
        mock_sig.id = uuid4()
        mock_sig.user_id = uuid4()
        mock_sig.action = "approve"
        mock_sig.entity_type = "compound"
        mock_sig.entity_id = uuid4()
        mock_sig.signature_hash = "hash123"
        mock_sig.meaning = "I approve"
        mock_sig.timestamp = datetime(2024, 1, 1)
        mock_sig.verified_at = None
        
        mock_get_sigs.return_value = [mock_sig]
        
        entity_id = uuid4()
        response = client.get(f"/api/v1/signatures/compound/{entity_id}")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        assert data[0]["action"] == "approve"

