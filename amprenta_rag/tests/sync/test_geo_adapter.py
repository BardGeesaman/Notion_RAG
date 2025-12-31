"""Tests for GEO sync adapter."""

import pytest
from unittest.mock import MagicMock, AsyncMock, patch
from datetime import datetime, timezone
from uuid import uuid4

from amprenta_rag.sync.adapters.geo import GEOSyncAdapter
from amprenta_rag.models.repository import StudyMetadata


class TestGEOSyncAdapter:
    """Test cases for GEO sync adapter."""

    def test_geo_adapter_source_name(self):
        """Test that adapter has correct source name."""
        adapter = GEOSyncAdapter()
        assert adapter.source == "geo"

    def test_geo_adapter_batch_size(self):
        """Test that adapter has appropriate batch size."""
        adapter = GEOSyncAdapter()
        assert adapter.BATCH_SIZE == 50
        assert adapter.BATCH_SIZE < 100  # Smaller than ChEMBL due to heavier metadata

    @pytest.mark.asyncio
    async def test_geo_adapter_fetch_records_basic(self):
        """Test basic fetch_records functionality without since date."""
        mock_repo = MagicMock()
        mock_repo._search_geo.return_value = ["GSE12345", "GSE67890"]
        
        # Mock metadata objects
        mock_metadata1 = MagicMock()
        mock_metadata1.title = "Study 1"
        mock_metadata1.summary = "Summary 1"
        mock_metadata1.organism = ["Homo sapiens"]
        mock_metadata1.platform = "GPL123"
        mock_metadata1.num_samples = 10
        mock_metadata1.omics_type = "transcriptomics"
        mock_metadata1.disease = "cancer"
        mock_metadata1.raw_metadata = {"key1": "value1"}
        
        mock_metadata2 = MagicMock()
        mock_metadata2.title = "Study 2"
        mock_metadata2.summary = "Summary 2"
        mock_metadata2.organism = ["Mus musculus"]
        mock_metadata2.platform = "GPL456"
        mock_metadata2.num_samples = 20
        mock_metadata2.omics_type = "transcriptomics"
        mock_metadata2.disease = "diabetes"
        mock_metadata2.raw_metadata = {"key2": "value2"}
        
        mock_repo.fetch_study_metadata.side_effect = [mock_metadata1, mock_metadata2]
        
        adapter = GEOSyncAdapter(geo_repo=mock_repo)
        records = []
        
        async for record in adapter.fetch_records(since=None):
            records.append(record)
        
        assert len(records) == 2
        
        # Verify first record
        assert records[0]["external_id"] == "GSE12345"
        assert records[0]["entity_type"] == "study"
        assert records[0]["data"]["study_id"] == "GSE12345"
        assert records[0]["data"]["title"] == "Study 1"
        assert records[0]["data"]["organism"] == ["human"]
        assert records[0]["data"]["repository"] == "GEO"
        
        # Verify search was called with correct query
        mock_repo._search_geo.assert_called_once_with(
            '"Homo sapiens"[Organism] AND gse[Entry Type]',
            max_results=50
        )

    @pytest.mark.asyncio
    async def test_geo_adapter_fetch_records_with_since_date(self):
        """Test fetch_records with since date for incremental sync."""
        mock_repo = MagicMock()
        mock_repo._search_geo.return_value = ["GSE12345"]
        
        mock_metadata = MagicMock()
        mock_metadata.title = "Study 1"
        mock_metadata.summary = "Summary 1"
        mock_metadata.organism = ["Homo sapiens"]
        mock_metadata.platform = "GPL123"
        mock_metadata.num_samples = 10
        mock_metadata.omics_type = "transcriptomics"
        mock_metadata.disease = "cancer"
        mock_metadata.raw_metadata = {"key1": "value1"}
        
        mock_repo.fetch_study_metadata.return_value = mock_metadata
        
        adapter = GEOSyncAdapter(geo_repo=mock_repo)
        since_date = datetime(2023, 1, 15, tzinfo=timezone.utc)
        
        records = []
        async for record in adapter.fetch_records(since=since_date):
            records.append(record)
        
        assert len(records) == 1
        
        # Verify MDAT filter was applied
        expected_query = '"Homo sapiens"[Organism] AND gse[Entry Type] AND 2023/01/15:3000[MDAT]'
        mock_repo._search_geo.assert_called_once_with(expected_query, max_results=50)

    def test_geo_adapter_mdat_filter_format(self):
        """Test that MDAT filter is formatted correctly."""
        adapter = GEOSyncAdapter()
        
        # Test various date formats
        test_dates = [
            (datetime(2023, 1, 1), "2023/01/01"),
            (datetime(2023, 12, 31), "2023/12/31"),
            (datetime(2024, 6, 15), "2024/06/15"),
        ]
        
        for test_date, expected_str in test_dates:
            formatted = test_date.strftime('%Y/%m/%d')
            assert formatted == expected_str

    def test_geo_adapter_checksum_computation(self):
        """Test checksum computation for change detection."""
        adapter = GEOSyncAdapter()
        
        # Test with sample record
        record = {
            "data": {
                "study_id": "GSE12345",
                "title": "Test Study",
                "summary": "Test Summary"
            }
        }
        
        checksum1 = adapter.compute_checksum(record)
        checksum2 = adapter.compute_checksum(record)
        
        # Same data should produce same checksum
        assert checksum1 == checksum2
        assert len(checksum1) == 32  # MD5 hex length
        
        # Different data should produce different checksum
        record["data"]["title"] = "Modified Title"
        checksum3 = adapter.compute_checksum(record)
        assert checksum3 != checksum1

    def test_geo_adapter_normalize_organism(self):
        """Test organism normalization functionality."""
        adapter = GEOSyncAdapter()
        
        # Test standard mappings
        test_cases = [
            (["Homo sapiens"], ["human"]),
            (["Mus musculus"], ["mouse"]),
            (["Rattus norvegicus"], ["rat"]),
            (["Danio rerio"], ["zebrafish"]),
            (["Drosophila melanogaster"], ["fruit fly"]),
            (["Caenorhabditis elegans"], ["c. elegans"]),
            (["Saccharomyces cerevisiae"], ["yeast"]),
            (["Unknown organism"], ["Unknown organism"]),  # No mapping
            (["Homo sapiens", "Mus musculus"], ["human", "mouse"]),  # Multiple
        ]
        
        for input_organisms, expected in test_cases:
            result = adapter._normalize_organism(input_organisms)
            assert result == expected
        
        # Test single string input
        result = adapter._normalize_organism("Homo sapiens")
        assert result == ["human"]
        
        # Test empty/None input
        result = adapter._normalize_organism([])
        assert result == []
        
        result = adapter._normalize_organism(None)
        assert result == []

    def test_geo_adapter_to_sync_record_format(self):
        """Test conversion of StudyMetadata to sync record format."""
        adapter = GEOSyncAdapter()
        
        # Create mock metadata
        mock_metadata = MagicMock()
        mock_metadata.title = "Test Study Title"
        mock_metadata.summary = "Test study summary"
        mock_metadata.organism = ["Homo sapiens", "Mus musculus"]
        mock_metadata.platform = "GPL570"
        mock_metadata.num_samples = 25
        mock_metadata.omics_type = "transcriptomics"
        mock_metadata.disease = "cancer"
        mock_metadata.raw_metadata = {"original": "data"}
        
        record = adapter._to_sync_record("GSE12345", mock_metadata)
        
        # Verify record structure
        assert record["external_id"] == "GSE12345"
        assert record["entity_type"] == "study"
        assert "checksum" in record
        assert "data" in record
        
        # Verify data fields
        data = record["data"]
        assert data["study_id"] == "GSE12345"
        assert data["title"] == "Test Study Title"
        assert data["summary"] == "Test study summary"
        assert data["organism"] == ["human", "mouse"]
        assert data["platform"] == "GPL570"
        assert data["num_samples"] == 25
        assert data["omics_type"] == "transcriptomics"
        assert data["disease"] == "cancer"
        assert data["repository"] == "GEO"

    @pytest.mark.asyncio
    async def test_geo_adapter_handles_fetch_errors(self):
        """Test that adapter handles fetch errors gracefully."""
        mock_repo = MagicMock()
        mock_repo._search_geo.return_value = ["GSE12345", "GSE67890", "GSE99999"]
        
        # First study succeeds, second fails, third succeeds
        mock_metadata1 = MagicMock()
        mock_metadata1.title = "Study 1"
        mock_metadata1.summary = "Summary 1"
        mock_metadata1.organism = ["Homo sapiens"]
        mock_metadata1.platform = "GPL123"
        mock_metadata1.num_samples = 10
        mock_metadata1.omics_type = "transcriptomics"
        mock_metadata1.disease = "cancer"
        mock_metadata1.raw_metadata = {"key1": "value1"}
        
        mock_metadata3 = MagicMock()
        mock_metadata3.title = "Study 3"
        mock_metadata3.summary = "Summary 3"
        mock_metadata3.organism = ["Homo sapiens"]
        mock_metadata3.platform = "GPL789"
        mock_metadata3.num_samples = 15
        mock_metadata3.omics_type = "transcriptomics"
        mock_metadata3.disease = "diabetes"
        mock_metadata3.raw_metadata = {"key3": "value3"}
        
        mock_repo.fetch_study_metadata.side_effect = [
            mock_metadata1,
            Exception("Network error"),
            mock_metadata3
        ]
        
        adapter = GEOSyncAdapter(geo_repo=mock_repo)
        records = []
        
        async for record in adapter.fetch_records(since=None):
            records.append(record)
        
        # Should get 2 records (first and third), skipping the failed one
        assert len(records) == 2
        assert records[0]["external_id"] == "GSE12345"
        assert records[1]["external_id"] == "GSE99999"

    def test_geo_adapter_composes_with_repo(self):
        """Test that adapter properly composes with GEORepository."""
        # Test with custom repo
        custom_repo = MagicMock()
        adapter = GEOSyncAdapter(geo_repo=custom_repo)
        assert adapter._repo is custom_repo
        
        # Test with default repo
        with patch('amprenta_rag.sync.adapters.geo.GEORepository') as mock_geo_repo:
            adapter = GEOSyncAdapter()
            mock_geo_repo.assert_called_once()

    def test_geo_adapter_map_to_entity(self):
        """Test mapping of external record to local entity."""
        adapter = GEOSyncAdapter()
        
        record = {
            "external_id": "GSE12345",
            "data": {"study_id": "GSE12345", "title": "Test Study"}
        }
        
        mock_db_session = MagicMock()
        entity_type, entity_id = adapter.map_to_entity(record, mock_db_session)
        
        assert entity_type == "study"
        assert entity_id is None  # New entity
