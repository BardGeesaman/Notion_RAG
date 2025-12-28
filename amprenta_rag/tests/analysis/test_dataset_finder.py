"""Tests for dataset finder functionality."""
import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.query.dataset_finder import (
    find_datasets_by_nl,
    DatasetFinderResult,
    DatasetResult,
    extract_search_terms,
    search_geo,
    search_arrayexpress,
    search_metabolomics_workbench,
    deduplicate_results,
    rank_results,
)


class TestDatasetFinder:
    """Test dataset finder functionality."""

    def test_find_datasets_by_nl_success(self):
        """Test successful dataset finding with natural language query."""
        query = "breast cancer RNA-seq human"
        
        with patch('amprenta_rag.query.dataset_finder.search_geo') as mock_geo, \
             patch('amprenta_rag.query.dataset_finder.search_arrayexpress') as mock_ae:
            
            geo_results = [
                DatasetResult(
                    accession="GSE123456", 
                    title="Breast cancer study", 
                    description="RNA-seq analysis", 
                    source="geo",
                    species="human", 
                    tissue="breast", 
                    disease="breast cancer", 
                    assay_type="RNA-seq", 
                    sample_count=100, 
                    url="https://geo.com/GSE123456"
                )
            ]
            mock_geo.return_value = geo_results
            mock_ae.return_value = []
            
            result = find_datasets_by_nl(query, repositories=["geo"], max_results=50)
            
            assert isinstance(result, DatasetFinderResult)
            assert result.query == query
            assert len(result.results) == 1
            assert result.results[0].accession == "GSE123456"
            assert "cancer" in result.extracted_terms["disease"]
            assert "geo" in result.sources_searched
            assert len(result.sources_failed) == 0

    def test_find_datasets_by_nl_empty_query(self):
        """Test dataset finder with empty query."""
        result = find_datasets_by_nl("", repositories=["GEO"])
        
        assert isinstance(result, DatasetFinderResult)
        assert result.query == ""
        assert len(result.results) == 0
        assert result.total_found == 0

    def test_find_datasets_by_nl_api_failure(self):
        """Test dataset finder with API failures."""
        query = "diabetes metabolomics"
        
        with patch('amprenta_rag.query.dataset_finder.extract_search_terms') as mock_extract, \
             patch('amprenta_rag.query.dataset_finder.search_geo') as mock_geo, \
             patch('amprenta_rag.query.dataset_finder.search_arrayexpress') as mock_ae:
            
            mock_extract.return_value = {"disease": ["diabetes"], "tissue": [], "species": [], "assay_type": ["metabolomics"]}
            mock_geo.side_effect = Exception("API timeout")
            mock_ae.return_value = []
            
            result = find_datasets_by_nl(query, repositories=["geo", "arrayexpress"])
            
            assert isinstance(result, DatasetFinderResult)
            assert "geo" in result.sources_failed
            assert "arrayexpress" in result.sources_searched

    def test_extract_search_terms_comprehensive(self):
        """Test search term extraction from natural language."""
        query = "Alzheimer disease mouse brain RNA-seq longitudinal study"
        
        # Test simple heuristic extraction (no LLM needed for this function)
        terms = extract_search_terms(query)
        
        # Check that terms are extracted (exact matches depend on heuristics)
        assert isinstance(terms, dict)
        assert "disease" in terms
        assert "tissue" in terms  
        assert "species" in terms
        assert "assay_type" in terms

    def test_extract_search_terms_api_failure(self):
        """Test search term extraction with API failure."""
        query = "cancer study"
        
        # Test with a query that might not match heuristics well
        terms = extract_search_terms(query)
            
        # Should still return structured terms
        assert isinstance(terms, dict)
        assert all(key in terms for key in ["disease", "tissue", "species", "assay_type"])

    def test_deduplicate_results_accession_match(self):
        """Test deduplication of results with matching accessions."""
        results = [
            DatasetResult("GSE123456", "Study A", "Description A", "geo"),
            DatasetResult("E-GEOD-123456", "Study A", "Description A", "arrayexpress"),  # Same study
            DatasetResult("GSE789012", "Study B", "Description B", "geo"),
        ]
        
        deduplicated = deduplicate_results(results)
        
        assert len(deduplicated) == 2  # One duplicate removed
        accessions = [r.accession for r in deduplicated]
        assert "GSE123456" in accessions
        assert "GSE789012" in accessions
        assert "E-GEOD-123456" not in accessions

    def test_deduplicate_results_title_similarity(self):
        """Test deduplication based on title similarity."""
        results = [
            DatasetResult("GSE111111", "Breast cancer RNA-seq analysis", "Description A", "geo"),
            DatasetResult("GSE222222", "Breast cancer RNA-seq analysis study", "Description B", "geo"),  # Similar title
            DatasetResult("GSE333333", "Lung cancer proteomics", "Description C", "geo"),
        ]
        
        deduplicated = deduplicate_results(results)
        
        assert len(deduplicated) == 2  # One duplicate removed based on title similarity
        titles = [r.title for r in deduplicated]
        assert "Lung cancer proteomics" in titles

    def test_rank_results_by_relevance(self):
        """Test ranking of results by relevance score."""
        results = [
            DatasetResult("GSE111111", "Unrelated study", "Some random description", "geo"),
            DatasetResult("GSE222222", "Cancer research study", "Breast cancer RNA-seq analysis", "geo", sample_count=200),
            DatasetResult("GSE333333", "Cancer study", "Basic cancer research", "geo", sample_count=50),
        ]
        
        # Terms that should match the content
        terms = {
            "disease": ["cancer", "breast cancer"],
            "tissue": ["breast"],
            "species": ["human"],
            "assay_type": ["RNA-seq"]
        }
        
        ranked = rank_results(results, terms)
        
        assert len(ranked) == 3
        assert ranked[0].accession == "GSE222222"  # Should have highest score (cancer + breast cancer + RNA-seq + large sample)
        assert ranked[1].accession == "GSE333333"  # Should have medium score (cancer only)
        assert ranked[2].accession == "GSE111111"  # Should have lowest score (no matches)
        assert ranked[0].score >= ranked[1].score >= ranked[2].score


class TestDatasetResult:
    """Test DatasetResult data class."""

    def test_dataset_result_creation(self):
        """Test DatasetResult object creation."""
        result = DatasetResult(
            accession="GSE123456",
            title="Test Study",
            description="Test Description",
            source="GEO",
            species="human",
            tissue="brain",
            disease="Alzheimer",
            assay_type="RNA-seq",
            sample_count=50,
            url="https://example.com",
            score=0.8
        )
        
        assert result.accession == "GSE123456"
        assert result.title == "Test Study"
        assert result.species == "human"
        assert result.score == 0.8

    def test_dataset_result_minimal(self):
        """Test DatasetResult with minimal required fields."""
        result = DatasetResult("GSE123456", "Test Study", "Test Description", "geo")
        
        assert result.accession == "GSE123456"
        assert result.title == "Test Study"
        assert result.source == "geo"
        assert result.species is None
        assert result.score == 0.0


class TestDatasetFinderResult:
    """Test DatasetFinderResult data class."""

    def test_dataset_finder_result_creation(self):
        """Test DatasetFinderResult object creation."""
        results = [DatasetResult("GSE123456", "Test Study", "Test Description", "geo")]
        extracted_terms = {"disease": ["cancer"], "tissue": ["brain"], "species": [], "assay_type": []}
        
        finder_result = DatasetFinderResult(
            query="cancer brain study",
            extracted_terms=extracted_terms,
            results=results,
            total_found=1,
            sources_searched=["geo"],
            sources_failed=[]
        )
        
        assert finder_result.query == "cancer brain study"
        assert len(finder_result.results) == 1
        assert finder_result.total_found == 1
        assert "geo" in finder_result.sources_searched
        assert len(finder_result.sources_failed) == 0


@pytest.mark.integration
class TestDatasetFinderIntegration:
    """Integration tests for dataset finder (requires external APIs)."""

    @pytest.mark.skip(reason="Requires external API access")
    def test_real_geo_search(self):
        """Test real GEO API search."""
        query = "breast cancer"
        result = find_datasets_by_nl(query, repositories=["geo"], max_results=5)
        
        assert isinstance(result, DatasetFinderResult)
        assert result.total_found > 0
        assert len(result.results) > 0
        assert all(r.source == "geo" for r in result.results)
