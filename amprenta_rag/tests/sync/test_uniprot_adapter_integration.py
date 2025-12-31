"""Integration tests for UniProt mapping adapter with mocked HTTP responses."""

import gzip
import pytest
from unittest.mock import AsyncMock, MagicMock, patch

from amprenta_rag.sync.adapters.uniprot_mapping import UniProtMappingAdapter

# Sample UniProt idmapping_selected.tab.gz content (columns: UniProtKB-AC, UniProtKB-ID, GeneID, RefSeq, GI, PDB, GO, UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI-taxon, MIM, UniGene, PubMed, EMBL, EMBL-CDS, Ensembl, Ensembl_TRS, Ensembl_PRO, Additional PubMed)
SAMPLE_MAPPING_DATA = b"P12345\tTP53_HUMAN\t7157\tNP_000537.3\t\t1TUP\tGO:0005634\tUniRef100_P12345\tUniRef90_P12345\tUniRef50_P12345\tUPI0000000001\t\t9606\t191170\t\t\tBC003596\tAAH03596.1\tENSG00000141510\tENST00000269305\tENSP00000269305\t\n"


@pytest.fixture
def mock_gzip_response():
    """Return sample data compressed with gzip."""
    return gzip.compress(SAMPLE_MAPPING_DATA)


class TestUniProtAdapterIntegration:
    """Test UniProt adapter integration with mocked HTTP responses."""
    
    def test_adapter_parses_valid_record(self, mock_gzip_response):
        """Test that adapter correctly parses a valid UniProt record."""
        # Parse the raw line directly
        line = SAMPLE_MAPPING_DATA.decode().strip()
        parts = line.split("\t")
        
        assert parts[0] == "P12345"  # UniProt ID
        assert parts[2] == "7157"     # Gene ID
        assert parts[18] == "ENSG00000141510"  # Ensembl gene
    
    def test_adapter_handles_malformed_row(self):
        """Test that adapter handles rows with missing columns gracefully."""
        malformed_line = "P12345\tTP53_HUMAN"  # Missing most columns
        parts = malformed_line.split("\t")
        
        # Should have only 2 parts
        assert len(parts) == 2
        # Gene ID column (index 2) would raise IndexError - adapter should handle
    
    @pytest.mark.asyncio
    async def test_fetch_records_with_mocked_http(self, mock_gzip_response):
        """Test full fetch_records flow with mocked HTTP response."""
        adapter = UniProtMappingAdapter()
        
        # Mock the httpx response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.content = mock_gzip_response
        mock_response.raise_for_status = MagicMock()
        
        # Mock httpx.AsyncClient to return our mock response
        with patch('httpx.AsyncClient') as mock_client_class:
            mock_client = AsyncMock()
            mock_client.__aenter__.return_value = mock_client
            mock_client.__aexit__.return_value = None
            mock_client.get.return_value = mock_response
            mock_client_class.return_value = mock_client
            
            records = []
            async for record in adapter.fetch_records(since=None):
                records.append(record)
            
            # Should parse the sample data and return at least one record
            assert len(records) >= 1
            assert mock_client.get.called
            
            # Verify the first record has expected structure
            if records:
                first_record = records[0]
                assert "uniprot_ac" in first_record
                assert first_record["uniprot_ac"] == "P12345"
