"""Tests for metadata enrichment functionality."""
import pytest
from unittest.mock import MagicMock, patch
from uuid import uuid4

from amprenta_rag.services.metadata_enrichment import (
    enrich_dataset_metadata,
    EnrichmentResult,
    _prepare_text_content,
    _merge_metadata_with_dataset,
    batch_enrich_datasets,
    get_enrichment_status,
)


class TestMetadataEnrichment:
    """Test metadata enrichment functionality."""

    def test_enrich_dataset_metadata_success(self):
        """Test successful dataset metadata enrichment."""
        dataset_id = uuid4()
        
        with patch('amprenta_rag.services.metadata_enrichment.db_session') as mock_db_session:
            
            # Mock database objects
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            
            mock_dataset = MagicMock()
            mock_dataset.id = dataset_id
            mock_dataset.title = "Test Dataset"
            mock_dataset.description = "A test dataset for enrichment"
            mock_dataset.metadata = {"existing_field": "value"}
            mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
            
            # Mock the enrichment dependencies
            with patch('amprenta_rag.services.metadata_enrichment._prepare_text_content') as mock_prepare, \
                 patch('amprenta_rag.services.metadata_enrichment.extract_semantic_metadata_with_llm') as mock_llm, \
                 patch('amprenta_rag.services.metadata_enrichment._merge_metadata_with_dataset') as mock_merge:
                
                mock_prepare.return_value = "Test dataset description for LLM processing"
                mock_llm.return_value = {"sample_count": 100, "study_design": "case-control"}
                mock_merge.return_value = ["sample_count", "study_design"]
                
                result = enrich_dataset_metadata(dataset_id)
                
                assert isinstance(result, EnrichmentResult)
                assert result.dataset_id == dataset_id
                assert result.success is True
                assert "sample_count" in result.enriched_fields
                assert result.extracted_metadata["sample_count"] == 100
                assert result.error_message is None

    def test_enrich_dataset_metadata_not_found(self):
        """Test metadata enrichment with dataset not found."""
        dataset_id = uuid4()
        
        with patch('amprenta_rag.services.metadata_enrichment.db_session') as mock_db_session:
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            mock_db.query.return_value.filter.return_value.first.return_value = None
            
            result = enrich_dataset_metadata(dataset_id)
            
            assert isinstance(result, EnrichmentResult)
            assert result.success is False
            assert "not found" in result.error_message.lower()
            assert len(result.enriched_fields) == 0

    def test_enrich_dataset_metadata_extraction_failure(self):
        """Test metadata enrichment with LLM extraction failure."""
        dataset_id = uuid4()
        
        with patch('amprenta_rag.services.metadata_enrichment.db_session') as mock_db_session, \
             patch('amprenta_rag.services.metadata_enrichment._prepare_text_content') as mock_prepare, \
             patch('amprenta_rag.services.metadata_enrichment.extract_semantic_metadata_with_llm') as mock_llm:
            
            # Mock database objects
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            
            mock_dataset = MagicMock()
            mock_dataset.id = dataset_id
            mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
            
            # Mock the enrichment pipeline with LLM failure
            mock_prepare.return_value = "Test dataset description"
            mock_llm.side_effect = Exception("LLM API error")
            
            result = enrich_dataset_metadata(dataset_id)
            
            assert isinstance(result, EnrichmentResult)
            assert result.success is False
            assert "LLM API error" in result.error_message

    def test_batch_enrich_datasets(self):
        """Test batch dataset enrichment."""
        dataset_ids = [uuid4(), uuid4()]
        
        with patch('amprenta_rag.services.metadata_enrichment.enrich_dataset_metadata') as mock_enrich:
            mock_results = [
                EnrichmentResult(
                    dataset_id=dataset_ids[0],
                    success=True,
                    enriched_fields=["sample_count"],
                    extracted_metadata={"sample_count": 100},
                    error_message=None,
                    processing_time_seconds=2.0
                ),
                EnrichmentResult(
                    dataset_id=dataset_ids[1],
                    success=True,
                    enriched_fields=["study_design"],
                    extracted_metadata={"study_design": "cohort"},
                    error_message=None,
                    processing_time_seconds=1.5
                )
            ]
            mock_enrich.side_effect = mock_results
            
            results = batch_enrich_datasets(dataset_ids)
            
            assert len(results) == 2
            assert all(r.success for r in results)
            assert results[0].extracted_metadata["sample_count"] == 100

    def test_get_enrichment_status(self):
        """Test getting enrichment status for a dataset."""
        dataset_id = uuid4()
        
        with patch('amprenta_rag.services.metadata_enrichment.db_session') as mock_db_session:
            mock_db = MagicMock()
            mock_db_session.return_value.__enter__.return_value = mock_db
            
            mock_dataset = MagicMock()
            mock_dataset.id = dataset_id
            mock_dataset.metadata = {
                "enrichment_timestamp": 1640995200.0,
                "enrichment_version": "1.0",
                "sample_count": 100,
                "study_design": "cohort"
            }
            mock_db.query.return_value.filter.return_value.first.return_value = mock_dataset
            
            status = get_enrichment_status(dataset_id)
            
            assert status["enriched"] is True
            assert len(status["available_fields"]) == 4


class TestEnrichmentResult:
    """Test EnrichmentResult data class."""

    def test_enrichment_result_success(self):
        """Test successful enrichment result creation."""
        dataset_id = uuid4()
        extracted_metadata = {"sample_count": 100, "study_design": "cohort"}
        
        result = EnrichmentResult(
            dataset_id=dataset_id,
            success=True,
            enriched_fields=["sample_count", "study_design"],
            extracted_metadata=extracted_metadata,
            error_message=None,
            processing_time_seconds=2.5
        )
        
        assert result.dataset_id == dataset_id
        assert result.success is True
        assert len(result.enriched_fields) == 2
        assert result.extracted_metadata["sample_count"] == 100
        assert result.error_message is None
        assert result.processing_time_seconds == 2.5

    def test_enrichment_result_failure(self):
        """Test failed enrichment result creation."""
        dataset_id = uuid4()
        
        result = EnrichmentResult(
            dataset_id=dataset_id,
            success=False,
            enriched_fields=[],
            extracted_metadata={},
            error_message="Dataset not found",
            processing_time_seconds=0.1
        )
        
        assert result.dataset_id == dataset_id
        assert result.success is False
        assert len(result.enriched_fields) == 0
        assert result.error_message == "Dataset not found"


@pytest.mark.integration
class TestMetadataEnrichmentIntegration:
    """Integration tests for metadata enrichment (requires LLM API)."""

    @pytest.mark.skip(reason="Requires OpenAI API access")
    def test_real_llm_extraction(self):
        """Test real LLM metadata extraction."""
        dataset_text = """
        This study analyzes gene expression in 50 breast cancer patients
        using RNA-sequencing. Samples were collected at diagnosis and
        after 3 months of treatment. The study uses a longitudinal design.
        """
        
        metadata = _extract_enhanced_metadata(dataset_text)
        
        assert isinstance(metadata, dict)
        assert "sample_count" in metadata
        assert metadata["sample_count"] > 0
        assert "study_design" in metadata
