"""Tests for UniProt incremental sync."""

import pytest
from unittest.mock import patch, MagicMock, AsyncMock

from amprenta_rag.sync.adapters.uniprot_mapping import UniProtMappingAdapter


class TestUniProtIncremental:
    """Test UniProt incremental sync features."""
    
    def test_adapter_source(self):
        """Test adapter source name."""
        adapter = UniProtMappingAdapter()
        assert adapter.source == "uniprot_mapping"
    
    def test_get_sync_metadata_empty(self):
        """Test sync metadata when no previous sync."""
        adapter = UniProtMappingAdapter()
        metadata = adapter.get_sync_metadata()
        assert metadata["etag"] is None
        assert metadata["last_modified"] is None
    
    @pytest.mark.asyncio
    async def test_304_not_modified_handling(self):
        """Test 304 Not Modified response handling."""
        adapter = UniProtMappingAdapter()
        
        # Mock the HTTP client to return 304
        mock_response = MagicMock()
        mock_response.status_code = 304
        
        with patch("amprenta_rag.sync.adapters.uniprot_mapping.db_session") as mock_db:
            mock_session = MagicMock()
            mock_session.__enter__ = MagicMock(return_value=mock_session)
            mock_session.__exit__ = MagicMock(return_value=False)
            mock_log = MagicMock()
            mock_log.metadata_ = {"etag": "abc123", "last_modified": "Mon, 01 Jan 2025 00:00:00 GMT"}
            mock_session.query.return_value.filter.return_value.order_by.return_value.first.return_value = mock_log
            mock_db.return_value = mock_session
            
            with patch("httpx.AsyncClient") as mock_client:
                mock_client_instance = AsyncMock()
                mock_client_instance.get = AsyncMock(return_value=mock_response)
                mock_client.return_value.__aenter__ = AsyncMock(return_value=mock_client_instance)
                mock_client.return_value.__aexit__ = AsyncMock(return_value=None)
                
                records = [r async for r in adapter.fetch_records()]
                assert len(records) == 0  # Empty iterator for 304


class TestUniProtMappingTransform:
    """Test UniProt mapping transformation."""
    
    @pytest.mark.asyncio
    async def test_transform_record_with_gene_id(self):
        """Test transformation with gene ID present."""
        adapter = UniProtMappingAdapter()
        
        record = {
            "uniprot_ac": "P12345",
            "gene_id": "7157",
            "ensembl": "ENSG00000141510",
        }
        
        mappings = await adapter.transform_record(record)
        
        # Should create bidirectional mappings
        assert len(mappings) >= 2
        
        # Check gene -> uniprot mapping
        gene_mapping = next((m for m in mappings if m["source_type"] == "gene"), None)
        assert gene_mapping is not None
        assert gene_mapping["target_id"] == "P12345"
    
    @pytest.mark.asyncio
    async def test_transform_record_empty(self):
        """Test transformation with empty record."""
        adapter = UniProtMappingAdapter()
        
        record = {}
        mappings = await adapter.transform_record(record)
        
        assert len(mappings) == 0
    
    @pytest.mark.asyncio
    async def test_transform_record_no_uniprot_ac(self):
        """Test transformation without UniProt AC."""
        adapter = UniProtMappingAdapter()
        
        record = {"gene_id": "7157"}
        mappings = await adapter.transform_record(record)
        
        assert len(mappings) == 0


class TestUniProtSyncBatch:
    """Test UniProt batch syncing."""
    
    @pytest.mark.asyncio
    async def test_sync_batch_empty(self):
        """Test sync_batch with empty records."""
        adapter = UniProtMappingAdapter()
        
        with patch.object(adapter, 'save_records', return_value=0):
            result = await adapter.sync_batch([])
            assert result == 0
    
    @pytest.mark.asyncio
    async def test_sync_batch_transforms_all_records(self):
        """Test sync_batch transforms and saves all records."""
        adapter = UniProtMappingAdapter()
        
        records = [
            {"uniprot_ac": "P12345", "gene_id": "7157"},
            {"uniprot_ac": "Q54321", "gene_id": "672"},
        ]
        
        with patch.object(adapter, 'save_records', return_value=8) as mock_save:
            result = await adapter.sync_batch(records)
            
            assert result == 8
            mock_save.assert_called_once()
            # Verify mappings were passed (should be multiple per record)
            saved_mappings = mock_save.call_args[0][0]
            assert len(saved_mappings) >= 4  # At least 2 mappings per record
