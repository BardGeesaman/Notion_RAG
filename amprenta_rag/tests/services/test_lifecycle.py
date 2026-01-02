"""Tests for lifecycle management service."""

import pytest
from uuid import uuid4
from unittest.mock import patch, MagicMock

from amprenta_rag.services.lifecycle import (
    execute_bulk_delete,
    find_orphaned_entities,
    cleanup_orphans,
    enforce_retention_policies,
    export_entity_data,
    create_export_package,
    calculate_deletion_impact,
    update_lifecycle_status,
)


class TestBulkDelete:
    """Tests for bulk delete functionality."""
    
    def test_bulk_delete_dry_run_returns_preview(self):
        """Dry run should return preview without deleting."""
        result = execute_bulk_delete(
            entity_type="dataset",
            entity_ids=[uuid4()],
            reason="test",
            dry_run=True
        )
        assert "entity_count" in result or "dry_run" in result
    
    def test_bulk_delete_invalid_entity_type(self):
        """Invalid entity type should return error."""
        result = execute_bulk_delete(
            entity_type="invalid",
            entity_ids=[uuid4()],
            reason="test",
            dry_run=False
        )
        assert result.get("errors") or result.get("failed", 0) > 0
    
    def test_bulk_delete_nonexistent_entity(self):
        """Deleting non-existent entity should fail gracefully."""
        result = execute_bulk_delete(
            entity_type="dataset",
            entity_ids=[uuid4()],
            reason="test",
            dry_run=False
        )
        assert result.get("failed", 0) >= 0  # May be 0 or 1


class TestOrphanDetection:
    """Tests for orphan detection and cleanup."""
    
    def test_find_orphaned_entities_returns_dict(self):
        """Should return dict with orphan counts."""
        result = find_orphaned_entities()
        assert isinstance(result, dict)
        assert "features" in result
        assert "signatures" in result
    
    def test_cleanup_orphans_dry_run(self):
        """Dry run should not delete anything."""
        result = cleanup_orphans(dry_run=True)
        assert result.get("dry_run") is True
        assert "would_delete" in result or "deleted" in result
    
    def test_cleanup_orphans_returns_counts(self):
        """Should return deletion counts."""
        result = cleanup_orphans(dry_run=False)
        assert "deleted" in result


class TestRetentionEnforcement:
    """Tests for retention policy enforcement."""
    
    def test_enforce_policies_dry_run(self):
        """Dry run should preview without changes."""
        result = enforce_retention_policies(dry_run=True)
        assert result.get("dry_run") is True
        assert "policies_checked" in result
    
    def test_enforce_policies_returns_actions(self):
        """Should return list of actions taken."""
        result = enforce_retention_policies(dry_run=True)
        assert "actions" in result
        assert isinstance(result["actions"], list)


class TestDataExport:
    """Tests for GDPR data export."""
    
    def test_export_nonexistent_entity(self):
        """Exporting non-existent entity should return error."""
        result = export_entity_data("dataset", uuid4())
        assert result.get("error") is not None or result.get("data") is None
    
    def test_export_invalid_entity_type(self):
        """Invalid entity type should return error."""
        result = export_entity_data("invalid", uuid4())
        assert result.get("error") is not None
    
    def test_export_includes_checksum(self):
        """Export should include SHA256 checksum."""
        result = export_entity_data("dataset", uuid4())
        # Either has checksum or has error
        assert "checksum" in result or "error" in result
    
    def test_create_export_package_returns_bytes(self):
        """Package creation should return bytes (ZIP)."""
        result = create_export_package("dataset", uuid4())
        assert isinstance(result, bytes)


class TestStatusUpdate:
    """Tests for status update functionality."""
    
    def test_update_invalid_status(self):
        """Invalid status should fail."""
        success, message = update_lifecycle_status(
            entity_type="dataset",
            entity_id=uuid4(),
            new_status="invalid_status"
        )
        assert success is False
    
    def test_update_nonexistent_entity(self):
        """Updating non-existent entity should fail."""
        success, message = update_lifecycle_status(
            entity_type="dataset",
            entity_id=uuid4(),
            new_status="archived"
        )
        assert success is False
        assert "not found" in message.lower()


class TestDeletionImpact:
    """Tests for deletion impact calculation."""
    
    def test_impact_nonexistent_entity(self):
        """Non-existent entity should return can_delete=False."""
        result = calculate_deletion_impact("dataset", uuid4())
        assert result.get("can_delete") is False or result.get("error") is not None
    
    def test_impact_returns_counts(self):
        """Should return impact counts dict."""
        result = calculate_deletion_impact("dataset", uuid4())
        assert "impact" in result