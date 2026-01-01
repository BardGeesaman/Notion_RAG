"""Unit tests for lifecycle management service."""

from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.services.lifecycle import (
    get_entity_by_id,
    calculate_deletion_impact,
    update_lifecycle_status,
    bulk_update_status,
    bulk_delete_preview,
    execute_bulk_archive,
    ENTITY_MODELS,
)
from amprenta_rag.database.models import LifecycleStatus


class TestGetEntityById:
    """Tests for get_entity_by_id function."""
    
    def test_returns_none_for_unknown_entity_type(self):
        """Unknown entity types return None."""
        result = get_entity_by_id("unknown", uuid4())
        assert result is None
    
    def test_valid_entity_types_in_mapping(self):
        """Verify all expected entity types are in mapping."""
        expected = {"dataset", "experiment", "compound", "signature"}
        assert set(ENTITY_MODELS.keys()) == expected


class TestCalculateDeletionImpact:
    """Tests for calculate_deletion_impact function."""
    
    @patch("amprenta_rag.services.lifecycle.db_session")
    def test_returns_error_for_nonexistent_entity(self, mock_db_session):
        """Non-existent entity returns error."""
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = None
        mock_db_session.return_value.__enter__.return_value = mock_session
        
        result = calculate_deletion_impact("dataset", uuid4())
        
        assert result["error"] == "Entity not found"
        assert result["can_delete"] is False
    
    def test_invalid_entity_type(self):
        """Invalid entity type returns empty impact."""
        result = calculate_deletion_impact("invalid_type", uuid4())
        assert "error" in result or result["can_delete"] is False


class TestUpdateLifecycleStatus:
    """Tests for update_lifecycle_status function."""
    
    def test_rejects_invalid_status(self):
        """Invalid status values are rejected."""
        success, message = update_lifecycle_status(
            "dataset", uuid4(), "nonexistent_status"
        )
        assert success is False
        assert "Invalid status" in message
    
    def test_rejects_unknown_entity_type(self):
        """Unknown entity types are rejected."""
        success, message = update_lifecycle_status(
            "unknown_type", uuid4(), "active"
        )
        assert success is False
        assert "Unknown entity type" in message
    
    @patch("amprenta_rag.services.lifecycle.db_session")
    def test_returns_unchanged_for_same_status(self, mock_db_session):
        """Same status returns unchanged message."""
        mock_entity = MagicMock()
        mock_entity.lifecycle_status = "active"
        
        mock_session = MagicMock()
        mock_session.query.return_value.filter.return_value.first.return_value = mock_entity
        mock_db_session.return_value.__enter__.return_value = mock_session
        
        success, message = update_lifecycle_status("dataset", uuid4(), "active")
        
        assert success is True
        assert "unchanged" in message.lower()


class TestBulkUpdateStatus:
    """Tests for bulk_update_status function."""
    
    @patch("amprenta_rag.services.lifecycle.update_lifecycle_status")
    def test_aggregates_results_correctly(self, mock_update):
        """Bulk update aggregates success/failure counts."""
        mock_update.side_effect = [
            (True, "Updated"),
            (False, "Entity not found"),
            (True, "Updated"),
        ]
        
        entity_ids = [uuid4(), uuid4(), uuid4()]
        results = bulk_update_status("dataset", entity_ids, "archived")
        
        assert results["total"] == 3
        assert results["success"] == 2
        assert results["failed"] == 1
        assert len(results["errors"]) == 1


class TestBulkDeletePreview:
    """Tests for bulk_delete_preview function."""
    
    @patch("amprenta_rag.services.lifecycle.calculate_deletion_impact")
    def test_aggregates_impact_counts(self, mock_impact):
        """Preview aggregates impact counts across entities."""
        mock_impact.side_effect = [
            {"impact": {"features": 10, "signatures": 5}, "blocking_references": []},
            {"impact": {"features": 20, "signatures": 3}, "blocking_references": []},
        ]
        
        preview = bulk_delete_preview("dataset", [uuid4(), uuid4()])
        
        assert preview["entity_count"] == 2
        assert preview["total_impact"]["features"] == 30
        assert preview["total_impact"]["signatures"] == 8
        assert preview["can_proceed"] is True
    
    @patch("amprenta_rag.services.lifecycle.calculate_deletion_impact")
    def test_blocking_references_prevent_proceed(self, mock_impact):
        """Blocking references set can_proceed to False."""
        mock_impact.return_value = {
            "impact": {"features": 10},
            "blocking_references": [{"type": "ml_model", "id": "123"}],
        }
        
        preview = bulk_delete_preview("dataset", [uuid4()])
        
        assert preview["can_proceed"] is False
        assert len(preview["blocking_entities"]) == 1


class TestExecuteBulkArchive:
    """Tests for execute_bulk_archive function."""
    
    @patch("amprenta_rag.services.lifecycle.bulk_update_status")
    def test_calls_bulk_update_with_archived_status(self, mock_bulk):
        """Execute bulk archive uses archived status."""
        mock_bulk.return_value = {"total": 2, "success": 2, "failed": 0, "errors": []}
        
        entity_ids = [uuid4(), uuid4()]
        execute_bulk_archive("compound", entity_ids, "Cleanup", uuid4())
        
        mock_bulk.assert_called_once()
        call_args = mock_bulk.call_args
        assert call_args[1]["new_status"] == LifecycleStatus.ARCHIVED.value
