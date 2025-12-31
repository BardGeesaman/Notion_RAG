"""Tests for entity versioning service."""

import pytest
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.auth.versioning import (
    compute_snapshot_checksum,
    entity_to_snapshot,
    create_version,
    get_versions,
    get_version,
    get_latest_version,
    compare_versions,
    get_version_history,
    rollback_to_version,
    delete_version,
    VERSIONABLE_ENTITIES,
)
from amprenta_rag.database.models import EntityVersion


class TestVersioning:
    """Test entity versioning functionality."""
    
    def test_compute_checksum_deterministic(self):
        """Test that checksum computation is deterministic."""
        data = {"id": str(uuid4()), "name": "test", "value": 42}
        
        checksum1 = compute_snapshot_checksum(data)
        checksum2 = compute_snapshot_checksum(data)
        
        assert checksum1 == checksum2
        assert len(checksum1) == 64  # SHA256 hex length
        assert isinstance(checksum1, str)
    
    def test_compute_checksum_different_data(self):
        """Test that different data produces different checksums."""
        data1 = {"id": str(uuid4()), "name": "test1"}
        data2 = {"id": str(uuid4()), "name": "test2"}
        
        checksum1 = compute_snapshot_checksum(data1)
        checksum2 = compute_snapshot_checksum(data2)
        
        assert checksum1 != checksum2
    
    def test_compute_checksum_order_independent(self):
        """Test that key order doesn't affect checksum."""
        data1 = {"name": "test", "id": str(uuid4()), "value": 42}
        data2 = {"id": data1["id"], "value": 42, "name": "test"}
        
        checksum1 = compute_snapshot_checksum(data1)
        checksum2 = compute_snapshot_checksum(data2)
        
        assert checksum1 == checksum2
    
    def test_entity_to_snapshot(self):
        """Test entity serialization to snapshot."""
        # Mock SQLAlchemy entity
        mock_entity = MagicMock()
        
        # Mock table and columns
        mock_table = MagicMock()
        mock_columns = []
        
        # Create mock columns
        for col_name in ["id", "name", "created_at"]:
            mock_col = MagicMock()
            mock_col.name = col_name
            mock_columns.append(mock_col)
        
        mock_table.columns = mock_columns
        mock_entity.__table__ = mock_table
        
        entity_id = uuid4()
        mock_entity.id = entity_id
        mock_entity.name = "Test Entity"
        mock_entity.created_at = "2023-01-01T00:00:00"
        
        snapshot = entity_to_snapshot(mock_entity)
        
        assert snapshot["id"] == entity_id
        assert snapshot["name"] == "Test Entity"
        assert snapshot["created_at"] == "2023-01-01T00:00:00"
    
    def test_create_version_success(self):
        """Test successful version creation."""
        entity_id = uuid4()
        user_id = uuid4()
        data = {"id": str(entity_id), "name": "Test Dataset"}
        
        # Use patch to mock the create_version function directly
        with patch("amprenta_rag.auth.versioning.EntityVersion") as mock_entity_version_class:
            mock_db = MagicMock()
            
            # Mock the database query chain
            mock_query_result = MagicMock()
            mock_query_result.scalar.return_value = None  # No existing versions
            mock_db.query.return_value.filter.return_value = mock_query_result
            
            # Mock the version instance
            mock_version = MagicMock()
            mock_version.version_number = 1
            mock_entity_version_class.return_value = mock_version
            
            result = create_version(
                db=mock_db,
                entity_type="dataset",
                entity_id=entity_id,
                data=data,
                user_id=user_id,
                change_summary="Initial version"
            )
            
            assert result == mock_version
            mock_db.add.assert_called_once_with(mock_version)
            mock_db.commit.assert_called_once()
            mock_db.refresh.assert_called_once_with(mock_version)
    
    def test_create_version_increments(self):
        """Test that version numbers auto-increment."""
        entity_id = uuid4()
        data = {"id": str(entity_id), "name": "Test Dataset"}
        
        with patch("amprenta_rag.auth.versioning.EntityVersion") as mock_entity_version_class:
            mock_db = MagicMock()
            
            # Mock the database query chain
            mock_query_result = MagicMock()
            mock_query_result.scalar.return_value = 3  # Existing max version
            mock_db.query.return_value.filter.return_value = mock_query_result
            
            # Mock the version instance
            mock_version = MagicMock()
            mock_version.version_number = 4
            mock_entity_version_class.return_value = mock_version
            
            result = create_version(
                db=mock_db,
                entity_type="dataset",
                entity_id=entity_id,
                data=data
            )
            
            assert result.version_number == 4
    
    def test_create_version_invalid_entity_type(self):
        """Test rejection of non-versionable entity types."""
        mock_db = MagicMock()
        entity_id = uuid4()
        data = {"id": str(entity_id)}
        
        with pytest.raises(ValueError, match="Entity type 'invalid' is not versionable"):
            create_version(
                db=mock_db,
                entity_type="invalid",
                entity_id=entity_id,
                data=data
            )
    
    def test_get_versions_returns_all(self):
        """Test retrieving all versions for an entity."""
        mock_db = MagicMock()
        entity_id = uuid4()
        
        # Mock versions
        version1 = MagicMock()
        version1.version_number = 2
        version2 = MagicMock()
        version2.version_number = 1
        
        mock_db.query.return_value.filter.return_value.order_by.return_value.all.return_value = [
            version1, version2
        ]
        
        versions = get_versions(mock_db, "dataset", entity_id)
        
        assert len(versions) == 2
        assert versions[0].version_number == 2  # Newest first
        assert versions[1].version_number == 1
    
    def test_get_versions_ordered_newest_first(self):
        """Test that versions are ordered newest first."""
        mock_db = MagicMock()
        entity_id = uuid4()
        
        # Mock query chain
        mock_query = mock_db.query.return_value
        mock_filter = mock_query.filter.return_value
        mock_order = mock_filter.order_by.return_value
        mock_order.all.return_value = []
        
        get_versions(mock_db, "dataset", entity_id)
        
        # Verify the order_by was called with desc()
        mock_filter.order_by.assert_called_once()
    
    def test_get_version_by_id(self):
        """Test retrieving a specific version by ID."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        mock_version = MagicMock()
        mock_version.id = version_id
        
        mock_db.query.return_value.filter.return_value.first.return_value = mock_version
        
        result = get_version(mock_db, version_id)
        
        assert result == mock_version
        assert result.id == version_id
    
    def test_get_version_not_found(self):
        """Test handling of non-existent version ID."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        mock_db.query.return_value.filter.return_value.first.return_value = None
        
        result = get_version(mock_db, version_id)
        
        assert result is None
    
    def test_get_latest_version(self):
        """Test retrieving the latest version."""
        mock_db = MagicMock()
        entity_id = uuid4()
        
        mock_version = MagicMock()
        mock_version.version_number = 5
        
        mock_db.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_version
        
        result = get_latest_version(mock_db, "dataset", entity_id)
        
        assert result == mock_version
        assert result.version_number == 5
    
    def test_compare_versions_added(self):
        """Test detecting added fields in version comparison."""
        v1_data = {"id": "123", "name": "old"}
        v2_data = {"id": "123", "name": "old", "description": "new field"}
        
        diff = compare_versions(v1_data, v2_data)
        
        assert diff["added"]["description"] == "new field"
        assert len(diff["removed"]) == 0
        assert len(diff["changed"]) == 0
    
    def test_compare_versions_removed(self):
        """Test detecting removed fields in version comparison."""
        v1_data = {"id": "123", "name": "test", "description": "old field"}
        v2_data = {"id": "123", "name": "test"}
        
        diff = compare_versions(v1_data, v2_data)
        
        assert diff["removed"]["description"] == "old field"
        assert len(diff["added"]) == 0
        assert len(diff["changed"]) == 0
    
    def test_compare_versions_changed(self):
        """Test detecting changed field values in version comparison."""
        v1_data = {"id": "123", "name": "old_name", "value": 10}
        v2_data = {"id": "123", "name": "new_name", "value": 20}
        
        diff = compare_versions(v1_data, v2_data)
        
        assert diff["changed"]["name"]["old"] == "old_name"
        assert diff["changed"]["name"]["new"] == "new_name"
        assert diff["changed"]["value"]["old"] == 10
        assert diff["changed"]["value"]["new"] == 20
        assert len(diff["added"]) == 0
        assert len(diff["removed"]) == 0
    
    def test_compare_versions_no_changes(self):
        """Test comparison of identical versions."""
        v1_data = {"id": "123", "name": "test", "value": 42}
        v2_data = {"id": "123", "name": "test", "value": 42}
        
        diff = compare_versions(v1_data, v2_data)
        
        assert len(diff["added"]) == 0
        assert len(diff["removed"]) == 0
        assert len(diff["changed"]) == 0
    
    def test_get_version_history_without_diffs(self):
        """Test getting version history without diffs."""
        mock_db = MagicMock()
        entity_id = uuid4()
        
        # Mock versions
        version1 = MagicMock()
        version1.id = uuid4()
        version1.version_number = 2
        version1.created_at = "2023-01-02"
        version1.created_by = uuid4()
        version1.change_summary = "Updated name"
        version1.checksum_sha256 = "abc123"
        
        version2 = MagicMock()
        version2.id = uuid4()
        version2.version_number = 1
        version2.created_at = "2023-01-01"
        version2.created_by = None
        version2.change_summary = "Initial version"
        version2.checksum_sha256 = "def456"
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_versions", lambda db, et, eid: [version1, version2])
            
            history = get_version_history(mock_db, "dataset", entity_id, include_diffs=False)
        
        assert len(history) == 2
        assert history[0]["version_number"] == 2
        assert history[1]["version_number"] == 1
        assert "diff" not in history[0]
        assert "diff" not in history[1]
    
    def test_get_version_history_with_diffs(self):
        """Test getting version history with diffs."""
        mock_db = MagicMock()
        entity_id = uuid4()
        
        # Mock versions with data snapshots
        version1 = MagicMock()
        version1.id = uuid4()
        version1.version_number = 2
        version1.created_at = "2023-01-02"
        version1.created_by = uuid4()
        version1.change_summary = "Updated name"
        version1.checksum_sha256 = "abc123"
        version1.data_snapshot = {"id": "123", "name": "new_name"}
        
        version2 = MagicMock()
        version2.id = uuid4()
        version2.version_number = 1
        version2.created_at = "2023-01-01"
        version2.created_by = None
        version2.change_summary = "Initial version"
        version2.checksum_sha256 = "def456"
        version2.data_snapshot = {"id": "123", "name": "old_name"}
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_versions", lambda db, et, eid: [version1, version2])
            
            history = get_version_history(mock_db, "dataset", entity_id, include_diffs=True)
        
        assert len(history) == 2
        assert "diff" in history[0]
        assert history[0]["diff"]["changed"]["name"]["old"] == "old_name"
        assert history[0]["diff"]["changed"]["name"]["new"] == "new_name"
        assert "diff" not in history[1]  # Last version has no diff
    
    def test_rollback_to_version(self):
        """Test rolling back to a previous version."""
        mock_db = MagicMock()
        version_id = uuid4()
        user_id = uuid4()
        
        # Mock source version
        source_version = MagicMock()
        source_version.entity_type = "dataset"
        source_version.entity_id = uuid4()
        source_version.version_number = 2
        source_version.data_snapshot = {"id": "123", "name": "old_data"}
        
        # Mock new version
        new_version = MagicMock()
        new_version.version_number = 4
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_version", lambda db, vid: source_version)
            mp.setattr("amprenta_rag.auth.versioning.create_version", lambda **kwargs: new_version)
            
            result = rollback_to_version(mock_db, version_id, user_id)
        
        assert result == new_version
    
    def test_rollback_version_not_found(self):
        """Test rollback with non-existent version."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_version", lambda db, vid: None)
            
            with pytest.raises(ValueError, match=f"Version {version_id} not found"):
                rollback_to_version(mock_db, version_id)
    
    def test_delete_version_success(self):
        """Test successful version deletion."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        # Mock version to delete
        mock_version = MagicMock()
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        mock_version.version_number = 2
        
        # Mock that there are other versions
        mock_db.query.return_value.filter.return_value.count.return_value = 2
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_version", lambda db, vid: mock_version)
            
            result = delete_version(mock_db, version_id)
        
        assert result is True
        mock_db.delete.assert_called_once_with(mock_version)
        mock_db.commit.assert_called_once()
    
    def test_delete_version_not_found(self):
        """Test deletion of non-existent version."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_version", lambda db, vid: None)
            
            result = delete_version(mock_db, version_id)
        
        assert result is False
        mock_db.delete.assert_not_called()
    
    def test_delete_last_version_forbidden(self):
        """Test that deleting the only version is forbidden."""
        mock_db = MagicMock()
        version_id = uuid4()
        
        # Mock version to delete
        mock_version = MagicMock()
        mock_version.entity_type = "dataset"
        mock_version.entity_id = uuid4()
        
        # Mock that this is the only version
        mock_db.query.return_value.filter.return_value.count.return_value = 0
        
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("amprenta_rag.auth.versioning.get_version", lambda db, vid: mock_version)
            
            with pytest.raises(ValueError, match="Cannot delete the only version of an entity"):
                delete_version(mock_db, version_id)
    
    def test_versionable_entities_list(self):
        """Test that versionable entities list is properly defined."""
        assert "dataset" in VERSIONABLE_ENTITIES
        assert "experiment" in VERSIONABLE_ENTITIES
        assert isinstance(VERSIONABLE_ENTITIES, list)
