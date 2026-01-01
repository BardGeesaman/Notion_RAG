"""Tests for auto-conflict resolution."""

from uuid import uuid4
from unittest.mock import patch, MagicMock

from amprenta_rag.sync.conflict_resolver import (
    ResolutionStrategy,
    resolve_conflict,
    auto_resolve_pending_conflicts,
    _apply_strategy,
    DEFAULT_STRATEGIES,
)


class TestResolutionStrategies:
    """Test conflict resolution strategies."""
    
    def test_prefer_external_strategy(self):
        """Test PREFER_EXTERNAL strategy."""
        local = {"value": "local_data"}
        external = {"value": "external_data"}
        
        result = _apply_strategy(ResolutionStrategy.PREFER_EXTERNAL, local, external)
        assert result == external
    
    def test_prefer_local_strategy(self):
        """Test PREFER_LOCAL strategy."""
        local = {"value": "local_data"}
        external = {"value": "external_data"}
        
        result = _apply_strategy(ResolutionStrategy.PREFER_LOCAL, local, external)
        assert result == local
    
    def test_prefer_newest_external_wins(self):
        """Test PREFER_NEWEST when external is newer."""
        local = {"value": "old", "updated_at": "2024-01-01"}
        external = {"value": "new", "updated_at": "2025-01-01"}
        
        result = _apply_strategy(ResolutionStrategy.PREFER_NEWEST, local, external)
        assert result == external
    
    def test_prefer_newest_local_wins(self):
        """Test PREFER_NEWEST when local is newer."""
        local = {"value": "new", "updated_at": "2025-01-01"}
        external = {"value": "old", "updated_at": "2024-01-01"}
        
        result = _apply_strategy(ResolutionStrategy.PREFER_NEWEST, local, external)
        assert result == local
    
    def test_default_strategies_defined(self):
        """Test default strategies are defined for common conflict types."""
        assert "schema_mismatch" in DEFAULT_STRATEGIES
        assert "value_conflict" in DEFAULT_STRATEGIES
        assert "missing_field" in DEFAULT_STRATEGIES
        assert "type_mismatch" in DEFAULT_STRATEGIES


class TestResolveConflict:
    """Test resolve_conflict function."""
    
    def test_resolve_not_found(self):
        """Test resolving non-existent conflict."""
        with patch("amprenta_rag.sync.conflict_resolver.db_session") as mock_db:
            mock_session = MagicMock()
            mock_session.__enter__ = MagicMock(return_value=mock_session)
            mock_session.__exit__ = MagicMock(return_value=False)
            mock_session.query.return_value.filter.return_value.first.return_value = None
            mock_db.return_value = mock_session
            
            result = resolve_conflict(uuid4())
            assert result is False
    
    def test_resolve_already_resolved(self):
        """Test resolving already-resolved conflict."""
        with patch("amprenta_rag.sync.conflict_resolver.db_session") as mock_db:
            mock_conflict = MagicMock()
            mock_conflict.resolution_status = "resolved"
            
            mock_session = MagicMock()
            mock_session.__enter__ = MagicMock(return_value=mock_session)
            mock_session.__exit__ = MagicMock(return_value=False)
            mock_session.query.return_value.filter.return_value.first.return_value = mock_conflict
            mock_db.return_value = mock_session
            
            result = resolve_conflict(uuid4())
            assert result is True


class TestAutoResolvePendingConflicts:
    """Test auto_resolve_pending_conflicts function."""
    
    def test_auto_resolve_empty(self):
        """Test auto-resolve with no pending conflicts."""
        with patch("amprenta_rag.sync.conflict_resolver.db_session") as mock_db:
            mock_session = MagicMock()
            mock_session.__enter__ = MagicMock(return_value=mock_session)
            mock_session.__exit__ = MagicMock(return_value=False)
            mock_session.query.return_value.filter.return_value.limit.return_value.all.return_value = []
            mock_db.return_value = mock_session
            
            result = auto_resolve_pending_conflicts()
            assert result == 0
