"""Tests for GEO sync adapter pagination improvements."""

import pytest
from datetime import datetime, timezone
from unittest.mock import patch

from amprenta_rag.sync.adapters.geo import GEOSyncAdapter


class TestGEOPagination:
    """Test GEO adapter pagination features."""
    
    def test_adapter_source(self):
        """Test adapter source name."""
        adapter = GEOSyncAdapter()
        assert adapter.source == "geo"
    
    def test_page_size_config(self):
        """Test page size configuration."""
        adapter = GEOSyncAdapter()
        assert adapter.PAGE_SIZE == 50
    
    def test_max_pages_config(self):
        """Test max pages configuration."""
        adapter = GEOSyncAdapter()
        assert adapter.MAX_PAGES == 20
    
    def test_rate_limit_config(self):
        """Test P1 rate limiting configuration."""
        adapter = GEOSyncAdapter()
        assert adapter.PAGE_RATE_LIMIT == 0.5
    
    def test_max_records_calculation(self):
        """Test maximum records per sync."""
        adapter = GEOSyncAdapter()
        max_records = adapter.PAGE_SIZE * adapter.MAX_PAGES
        assert max_records == 1000
    
    @pytest.mark.asyncio
    async def test_fetch_records_empty_results(self):
        """Test empty results handling."""
        adapter = GEOSyncAdapter()
        
        with patch.object(adapter._repo, '_search_geo', return_value=[]):
            with patch.object(adapter._repo, 'fetch_study_metadata'):
                records = [r async for r in adapter.fetch_records()]
                assert len(records) == 0
    
    @pytest.mark.asyncio
    async def test_incremental_query_format(self):
        """Test incremental query includes MDAT filter."""
        adapter = GEOSyncAdapter()
        since = datetime(2025, 1, 1, tzinfo=timezone.utc)
        
        with patch.object(adapter._repo, '_search_geo', return_value=[]) as mock_search:
            [r async for r in adapter.fetch_records(since=since)]
            
            # Verify MDAT filter was included in query
            call_args = mock_search.call_args
            query = call_args[0][0] if call_args[0] else call_args[1].get('query', '')
            assert "MDAT" in query or "2025/01/01" in query
    
    def test_compute_checksum(self):
        """Test checksum computation."""
        adapter = GEOSyncAdapter()
        record = {"data": {"study_id": "GSE12345", "title": "Test Study"}}
        checksum = adapter.compute_checksum(record)
        
        assert isinstance(checksum, str)
        assert len(checksum) == 32  # MD5 hex length
    
    def test_map_to_entity(self):
        """Test entity mapping returns study type."""
        adapter = GEOSyncAdapter()
        entity_type, entity_id = adapter.map_to_entity({}, None)
        
        assert entity_type == "study"
        assert entity_id is None
