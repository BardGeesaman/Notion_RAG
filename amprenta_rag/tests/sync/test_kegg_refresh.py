"""Tests for KEGG cache refresh adapter."""

import pytest
from unittest.mock import patch, MagicMock

from amprenta_rag.sync.adapters.kegg_refresh import KEGGRefreshAdapter


class TestKEGGRefreshAdapter:
    """Test KEGG cache refresh adapter."""
    
    def test_source_name(self):
        """Test adapter source name."""
        adapter = KEGGRefreshAdapter()
        assert adapter.source == "kegg_refresh"
    
    def test_refresh_window_days(self):
        """Test refresh window configuration."""
        adapter = KEGGRefreshAdapter()
        assert adapter.REFRESH_WINDOW_DAYS == 7
    
    def test_rate_limit_seconds(self):
        """Test rate limiting configuration."""
        adapter = KEGGRefreshAdapter()
        assert adapter.RATE_LIMIT_SECONDS == 0.35
    
    @pytest.mark.asyncio
    async def test_fetch_records_empty_when_no_expiring(self):
        """Test empty iterator when no expiring mappings."""
        adapter = KEGGRefreshAdapter()
        
        with patch("amprenta_rag.sync.adapters.kegg_refresh.db_session") as mock_db:
            mock_session = MagicMock()
            mock_session.__enter__ = MagicMock(return_value=mock_session)
            mock_session.__exit__ = MagicMock(return_value=False)
            mock_session.query.return_value.filter.return_value.limit.return_value.all.return_value = []
            mock_db.return_value = mock_session
            
            records = [r async for r in adapter.fetch_records()]
            assert len(records) == 0
    
    def test_compute_checksum(self):
        """Test checksum computation."""
        adapter = KEGGRefreshAdapter()
        
        record = {"source_id": "TP53", "new_target_id": "hsa:7157"}
        checksum = adapter.compute_checksum(record)
        
        assert isinstance(checksum, str)
        assert len(checksum) == 32  # MD5 hex length
    
    def test_map_to_entity(self):
        """Test entity mapping."""
        adapter = KEGGRefreshAdapter()
        
        entity_type, entity_id = adapter.map_to_entity({}, None)
        assert entity_type == "id_mapping"
        assert entity_id is None
