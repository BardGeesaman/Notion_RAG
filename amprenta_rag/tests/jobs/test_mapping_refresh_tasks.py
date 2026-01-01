"""Tests for mapping refresh Celery tasks."""

from unittest.mock import patch, MagicMock

from amprenta_rag.jobs.tasks.mapping_refresh import (
    refresh_kegg_cache_task,
    auto_resolve_sync_conflicts_task,
    cleanup_expired_mappings_task,
)


class TestRefreshKEGGCacheTask:
    """Test KEGG cache refresh task."""
    
    def test_task_name(self):
        """Test task has correct name."""
        assert refresh_kegg_cache_task.name == "refresh_kegg_cache"
    
    @patch("amprenta_rag.services.id_mapping_service.log_refresh_start")
    @patch("amprenta_rag.services.id_mapping_service.log_refresh_complete")
    @patch("amprenta_rag.sync.adapters.kegg_refresh.KEGGRefreshAdapter")
    def test_task_success(self, mock_adapter_cls, mock_log_complete, mock_log_start):
        """Test successful KEGG refresh."""
        # Setup mocks
        mock_log_entry = MagicMock()
        mock_log_entry.id = "test-log-id"
        mock_log_start.return_value = mock_log_entry
        
        mock_adapter = MagicMock()
        
        async def mock_fetch():
            yield {"refreshed": True}
            yield {"refreshed": True}
            yield {"refreshed": False}
        
        mock_adapter.fetch_records = mock_fetch
        mock_adapter_cls.return_value = mock_adapter
        
        # Execute task
        result = refresh_kegg_cache_task()
        
        # Verify
        assert result["status"] == "success"
        assert result["refreshed_count"] == 2
        mock_log_start.assert_called_once_with("kegg_refresh")
        mock_log_complete.assert_called_once()


class TestAutoResolveSyncConflictsTask:
    """Test auto-resolve sync conflicts task."""
    
    def test_task_name(self):
        """Test task has correct name."""
        assert auto_resolve_sync_conflicts_task.name == "auto_resolve_sync_conflicts"
    
    @patch("amprenta_rag.sync.conflict_resolver.auto_resolve_pending_conflicts")
    def test_task_success(self, mock_resolve):
        """Test successful conflict resolution."""
        mock_resolve.return_value = 5
        
        result = auto_resolve_sync_conflicts_task()
        
        assert result["status"] == "success"
        assert result["resolved_count"] == 5
        mock_resolve.assert_called_once_with(source=None, limit=100)
    
    @patch("amprenta_rag.sync.conflict_resolver.auto_resolve_pending_conflicts")
    def test_task_with_source_filter(self, mock_resolve):
        """Test conflict resolution with source filter."""
        mock_resolve.return_value = 3
        
        result = auto_resolve_sync_conflicts_task(source="kegg_refresh", limit=50)
        
        assert result["resolved_count"] == 3
        mock_resolve.assert_called_once_with(source="kegg_refresh", limit=50)


class TestCleanupExpiredMappingsTask:
    """Test cleanup expired mappings task."""
    
    def test_task_name(self):
        """Test task has correct name."""
        assert cleanup_expired_mappings_task.name == "cleanup_expired_mappings"
    
    @patch("amprenta_rag.services.id_mapping_service.cleanup_expired_mappings")
    def test_task_success(self, mock_cleanup):
        """Test successful cleanup."""
        mock_cleanup.return_value = 10
        
        result = cleanup_expired_mappings_task()
        
        assert result["status"] == "success"
        assert result["deleted_count"] == 10
        mock_cleanup.assert_called_once()
