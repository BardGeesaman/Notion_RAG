"""Tests for audit checksum extensions."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from amprenta_rag.auth.audit import (
    compute_checksum,
    get_audit_trail,
    log_data_change,
    verify_integrity,
)


class TestComputeChecksum:
    """Tests for compute_checksum."""

    def test_compute_checksum_deterministic(self):
        """Test that checksum is deterministic."""
        data = {"name": "test", "value": 123}
        
        checksum1 = compute_checksum(data)
        checksum2 = compute_checksum(data)
        
        assert checksum1 == checksum2
        assert len(checksum1) == 64  # SHA256 hex is 64 chars

    def test_compute_checksum_different_data(self):
        """Test that different data produces different checksums."""
        data1 = {"name": "test1"}
        data2 = {"name": "test2"}
        
        checksum1 = compute_checksum(data1)
        checksum2 = compute_checksum(data2)
        
        assert checksum1 != checksum2


class TestLogDataChange:
    """Tests for log_data_change."""

    def test_log_data_change_creates_record(self):
        """Test that log_data_change creates audit record."""
        mock_db = MagicMock()
        
        old_data = {"name": "old_value"}
        new_data = {"name": "new_value"}
        
        log_data_change(
            mock_db,
            entity_type="dataset",
            entity_id="test-id",
            action="update",
            actor_id=None,
            old_data=old_data,
            new_data=new_data,
        )
        
        assert mock_db.add.called
        assert mock_db.commit.called


class TestGetAuditTrail:
    """Tests for get_audit_trail."""

    def test_get_audit_trail_returns_list(self):
        """Test that get_audit_trail returns list."""
        mock_log = MagicMock()
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = [mock_log]
        
        result = get_audit_trail(mock_db, "dataset", "test-id")
        
        assert isinstance(result, list)
        assert len(result) == 1


class TestVerifyIntegrity:
    """Tests for verify_integrity."""

    def test_verify_integrity_no_audit(self):
        """Test integrity verification with no audit record."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = None
        
        result = verify_integrity(mock_db, "dataset", "test-id", {"name": "test"})
        
        assert result is True  # No audit record = pass

    def test_verify_integrity_matching_checksum(self):
        """Test integrity verification with matching checksum."""
        current_data = {"name": "test"}
        current_checksum = compute_checksum(current_data)
        
        mock_log = MagicMock()
        mock_log.new_checksum = current_checksum
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_log
        
        result = verify_integrity(mock_db, "dataset", "test-id", current_data)
        
        assert result is True

    def test_verify_integrity_mismatched_checksum(self):
        """Test integrity verification with mismatched checksum."""
        current_data = {"name": "modified"}
        
        mock_log = MagicMock()
        mock_log.new_checksum = "different_checksum_abc123"
        
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_log
        
        result = verify_integrity(mock_db, "dataset", "test-id", current_data)
        
        assert result is False

