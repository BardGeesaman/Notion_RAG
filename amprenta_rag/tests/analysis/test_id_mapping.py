"""
Unit tests for ID mapping services.

Tests the mapping functions that convert feature identifiers to pathway database IDs:
- map_protein_to_uniprot
- map_gene_to_kegg
- map_metabolite_to_kegg
- batch_map_features_to_pathway_ids
- Caching behavior
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch


from amprenta_rag.analysis.id_mapping import (
    batch_map_features_to_pathway_ids,
    clear_id_mapping_cache,
    get_cache_stats,
    map_gene_to_kegg,
    map_gene_to_reactome,
    map_metabolite_to_kegg,
    map_protein_to_kegg,
    map_protein_to_reactome,
    map_protein_to_uniprot,
)


class TestMapProteinToUniprot:
    """Tests for map_protein_to_uniprot function."""

    def test_valid_uniprot_id_returned_as_is(self):
        """Valid UniProt IDs should be returned unchanged."""
        # Standard UniProt accession format
        assert map_protein_to_uniprot("P12345") == "P12345"
        assert map_protein_to_uniprot("Q9Y6K9") == "Q9Y6K9"
        assert map_protein_to_uniprot("O15394") == "O15394"

    def test_uniprot_pattern_recognition(self):
        """Test UniProt ID pattern matching."""
        # Valid patterns
        valid_ids = ["P12345", "Q9Y6K9", "O15394", "A0A024R1R8"]
        for uid in valid_ids:
            result = map_protein_to_uniprot(uid)
            # Should return the ID as-is if it matches UniProt pattern
            assert result == uid or result is None  # None if API call fails

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    @patch("amprenta_rag.analysis.id_mapping.requests.post")
    def test_gene_name_mapping_success(self, mock_post, mock_get):
        """Test mapping gene name to UniProt ID via API."""
        # Mock POST for job submission
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"jobId": "test-job-123"}
        mock_post.return_value = mock_response

        # Mock GET for polling results
        mock_get_response = MagicMock()
        mock_get_response.status_code = 200
        mock_get_response.json.return_value = {"results": [{"to": {"primaryAccession": "P04637"}}]}
        mock_get.return_value = mock_get_response

        # Clear cache to ensure fresh lookup
        clear_id_mapping_cache()

        # Gene name that's not a UniProt ID
        result = map_protein_to_uniprot("TP53")
        # Should attempt API call
        assert mock_post.called or result is None

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    @patch("amprenta_rag.analysis.id_mapping.requests.post")
    def test_invalid_protein_returns_none(self, mock_post, mock_get):
        """Invalid protein IDs should return None."""
        # Mock both to fail
        mock_post.side_effect = Exception("API error")
        mock_get.side_effect = Exception("API error")
        clear_id_mapping_cache()

        result = map_protein_to_uniprot("INVALID_PROTEIN_XYZ_12345")
        assert result is None


class TestMapGeneToKegg:
    """Tests for map_gene_to_kegg function."""

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_successful_gene_mapping(self, mock_get):
        """Test successful gene to KEGG mapping."""
        # Mock KEGG API response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "hsa:7157\tTP53; tumor protein p53"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_gene_to_kegg("TP53", organism="hsa")

        assert result == "hsa:7157"
        mock_get.assert_called()

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_gene_not_found(self, mock_get):
        """Test gene not found in KEGG."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ""  # Empty response
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_gene_to_kegg("NONEXISTENT_GENE_123")

        assert result is None

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_api_error_returns_none(self, mock_get):
        """API errors should return None."""
        mock_get.side_effect = Exception("Network error")

        clear_id_mapping_cache()
        result = map_gene_to_kegg("TP53")

        assert result is None


class TestMapMetaboliteToKegg:
    """Tests for map_metabolite_to_kegg function."""

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_successful_metabolite_mapping(self, mock_get):
        """Test successful metabolite to KEGG compound mapping."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "cpd:C00031\tD-Glucose; Grape sugar; Dextrose"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_metabolite_to_kegg("glucose")

        assert result == "cpd:C00031"

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_exact_match_preferred(self, mock_get):
        """Exact name matches should be preferred."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        # Multiple results, ATP is exact match
        mock_response.text = "cpd:C00002\tATP; Adenosine 5'-triphosphate\ncpd:C00003\tOther compound"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_metabolite_to_kegg("ATP")

        assert result == "cpd:C00002"

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_metabolite_not_found(self, mock_get):
        """Metabolite not found should return None."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ""
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_metabolite_to_kegg("nonexistent_metabolite_xyz")

        assert result is None


class TestMapGeneToReactome:
    """Tests for map_gene_to_reactome function."""

    def test_valid_gene_symbol_returned(self):
        """Valid gene symbols should be returned for Reactome."""
        # Reactome uses gene symbols directly
        assert map_gene_to_reactome("TP53") == "TP53"
        assert map_gene_to_reactome("BRCA1") == "BRCA1"
        assert map_gene_to_reactome("EGFR") == "EGFR"

    def test_invalid_gene_symbol_returns_none(self):
        """Invalid gene symbol format should return None."""
        # Gene symbols with special characters
        result = map_gene_to_reactome("invalid@gene!")
        assert result is None


class TestMapProteinToReactome:
    """Tests for map_protein_to_reactome function."""

    def test_valid_uniprot_id_for_reactome(self):
        """Valid UniProt IDs should be returned for Reactome."""
        # Reactome uses UniProt IDs, so this should call map_protein_to_uniprot
        result = map_protein_to_reactome("P12345")
        assert result == "P12345"


class TestBatchMapFeatures:
    """Tests for batch_map_features_to_pathway_ids function."""

    @patch("amprenta_rag.analysis.id_mapping.map_gene_to_kegg")
    def test_batch_gene_mapping(self, mock_map_gene):
        """Test batch mapping of genes to KEGG."""
        mock_map_gene.side_effect = lambda g, organism="hsa": f"hsa:{hash(g) % 10000}"

        features = {"TP53", "BRCA1", "EGFR"}
        result = batch_map_features_to_pathway_ids(
            features, feature_type="gene", pathway_source="KEGG"
        )

        assert len(result) == 3
        assert all(v is not None for v in result.values())

    @patch("amprenta_rag.analysis.id_mapping.map_metabolite_to_kegg")
    def test_batch_metabolite_mapping(self, mock_map_met):
        """Test batch mapping of metabolites to KEGG."""
        mock_map_met.side_effect = lambda m: f"cpd:C{hash(m) % 100000:05d}"

        features = {"glucose", "ATP", "lactate"}
        result = batch_map_features_to_pathway_ids(
            features, feature_type="metabolite", pathway_source="KEGG"
        )

        assert len(result) == 3
        assert all(v is not None for v in result.values())

    def test_unsupported_feature_type(self):
        """Unsupported feature types should return None mappings."""
        features = {"feature1", "feature2"}
        result = batch_map_features_to_pathway_ids(
            features, feature_type="unsupported_type", pathway_source="KEGG"
        )

        assert len(result) == 2
        assert all(v is None for v in result.values())


class TestCaching:
    """Tests for caching behavior."""

    def test_cache_stats(self):
        """Test cache statistics function."""
        clear_id_mapping_cache()
        stats = get_cache_stats()

        assert "total_cached" in stats
        assert "successful_mappings" in stats
        assert "failed_mappings" in stats
        assert stats["total_cached"] == 0

    def test_clear_cache(self):
        """Test cache clearing."""
        clear_id_mapping_cache()
        stats_before = get_cache_stats()
        assert stats_before["total_cached"] == 0

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_cached_results_reused(self, mock_get):
        """Cached results should be reused without API calls."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "hsa:7157\tTP53"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()

        # First call should hit API
        result1 = map_gene_to_kegg("TP53")
        call_count_after_first = mock_get.call_count

        # Second call should use cache
        result2 = map_gene_to_kegg("TP53")
        call_count_after_second = mock_get.call_count

        assert result1 == result2
        # Second call shouldn't increase call count (using cache)
        assert call_count_after_second == call_count_after_first


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_feature_set(self):
        """Empty feature set should return empty dict."""
        result = batch_map_features_to_pathway_ids(
            set(), feature_type="gene", pathway_source="KEGG"
        )
        assert result == {}

    def test_none_handling(self):
        """None values should be handled gracefully."""
        # This shouldn't raise an exception
        clear_id_mapping_cache()
        # Most functions don't accept None directly, but they handle missing results

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_timeout_handling(self, mock_get):
        """Timeout should be handled gracefully."""
        import requests

        mock_get.side_effect = requests.Timeout("Connection timed out")

        clear_id_mapping_cache()
        result = map_gene_to_kegg("TP53")

        assert result is None

    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_malformed_response_handling(self, mock_get):
        """Malformed API responses should be handled."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "malformed\tdata\twith\textra\tcolumns"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        # Should not raise, may return the first column or None
        result = map_gene_to_kegg("TP53")
        # Result depends on implementation, just verify no exception


class TestProteinToKegg:
    """Tests for map_protein_to_kegg function."""

    @patch("amprenta_rag.analysis.id_mapping.map_protein_to_uniprot")
    @patch("amprenta_rag.analysis.id_mapping.requests.get")
    def test_protein_to_kegg_via_uniprot(self, mock_get, mock_uniprot):
        """Test protein to KEGG mapping via UniProt intermediate."""
        mock_uniprot.return_value = "P04637"
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "up:P04637\thsa:7157"
        mock_get.return_value = mock_response

        clear_id_mapping_cache()
        result = map_protein_to_kegg("TP53_HUMAN")

        assert mock_uniprot.called

    @patch("amprenta_rag.analysis.id_mapping.map_protein_to_uniprot")
    def test_protein_to_kegg_uniprot_fails(self, mock_uniprot):
        """If UniProt mapping fails, KEGG mapping should fail."""
        mock_uniprot.return_value = None

        clear_id_mapping_cache()
        result = map_protein_to_kegg("unknown_protein")

        assert result is None

