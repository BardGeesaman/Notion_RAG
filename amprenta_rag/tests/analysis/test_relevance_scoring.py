"""Tests for relevance and novelty scoring functionality."""
import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from uuid import uuid4

from amprenta_rag.analysis.relevance_scoring import (
    score_relevance,
    score_novelty,
    batch_score,
    RelevanceScore,
    NoveltyScore,
    ScoredItem,
    _get_item_embedding,
    _get_from_cache,
    _set_cache,
)


class TestRelevanceScoring:
    """Test relevance scoring functionality."""

    def test_score_relevance_success(self):
        """Test successful relevance scoring."""
        item = {
            "id": "dataset_123",
            "title": "Breast cancer RNA-seq study",
            "description": "Gene expression analysis in breast cancer patients",
            "species": "human",
            "assay_type": "RNA-seq"
        }
        
        context = {
            "diseases": ["breast cancer"],
            "targets": ["BRCA1", "BRCA2"],
            "species": ["human"],
            "assay_types": ["RNA-seq"]
        }
        
        criteria = {
            "disease_match": 0.4,
            "target_overlap": 0.3,
            "data_quality": 0.3
        }
        
        with patch('amprenta_rag.analysis.relevance_scoring.get_openai_client') as mock_client:
            mock_response = MagicMock()
            mock_response.choices[0].message.content = """{
                "overall_score": 0.85,
                "disease_match": 0.9,
                "target_overlap": 0.7,
                "data_quality": 0.8,
                "explanation": "High relevance due to exact disease match and compatible assay type"
            }"""
            mock_client.return_value.chat.completions.create.return_value = mock_response
            
            score = score_relevance(item, context, criteria)
            
            assert isinstance(score, RelevanceScore)
            assert score.item_id == "dataset_123"
            assert score.overall_score == 0.85
            assert score.disease_match == 0.9
            assert score.target_overlap == 0.7
            assert score.data_quality == 0.8
            assert "High relevance" in score.explanation

    def test_score_relevance_api_failure(self):
        """Test relevance scoring with LLM API failure."""
        item = {"id": "dataset_123", "title": "Test dataset"}
        context = {"diseases": ["cancer"]}
        
        with patch('amprenta_rag.analysis.relevance_scoring.get_openai_client') as mock_client:
            mock_client.return_value.chat.completions.create.side_effect = Exception("API error")
            
            score = score_relevance(item, context)
            
            assert isinstance(score, RelevanceScore)
            # On error, function returns 0.1 as default low score
            assert score.overall_score == 0.1
            assert "error" in score.explanation.lower()

    def test_score_novelty_success(self):
        """Test successful novelty scoring."""
        item = {
            "id": "new_dataset",
            "title": "Novel cancer biomarker study",
            "description": "Discovery of new biomarkers"
        }
        
        existing_items = [
            {
                "id": "existing_1",
                "title": "Traditional cancer study",
                "description": "Standard cancer analysis"
            },
            {
                "id": "existing_2", 
                "title": "Heart disease research",
                "description": "Cardiovascular analysis"
            }
        ]
        
        # Mock the internal dependencies, not score_novelty itself
        with patch('amprenta_rag.analysis.relevance_scoring._get_item_embedding') as mock_embed, \
             patch('amprenta_rag.analysis.relevance_scoring._get_novelty_explanation') as mock_explain:
            
            # Return different embeddings to simulate low similarity
            mock_embed.side_effect = [
                np.array([1.0, 0.0, 0.0]),  # item
                np.array([0.0, 1.0, 0.0]),  # existing_1
                np.array([0.0, 0.0, 1.0]),  # existing_2
            ]
            mock_explain.return_value = "Novel approach to biomarker discovery"
            
            score = score_novelty(item, existing_items)
            
            assert isinstance(score, NoveltyScore)
            assert score.item_id == "new_dataset"
            # With orthogonal embeddings, similarity should be 0, novelty should be 1.0
            assert score.novelty_score >= 0.9  # High novelty due to low similarity
            assert score.max_similarity <= 0.1  # Low similarity to existing items

    def test_score_novelty_high_similarity(self):
        """Test novelty scoring with high similarity to existing items."""
        item = {
            "id": "similar_dataset",
            "title": "Cancer study",
            "description": "Cancer analysis"
        }
        
        existing_items = [
            {
                "id": "existing_1",
                "title": "Cancer study protocol",
                "description": "Cancer analysis methods"
            }
        ]
        
        # Mock the internal dependencies to simulate high similarity
        with patch('amprenta_rag.analysis.relevance_scoring._get_item_embedding') as mock_embed, \
             patch('amprenta_rag.analysis.relevance_scoring._get_novelty_explanation') as mock_explain:
            
            # Return very similar embeddings to simulate high similarity
            mock_embed.side_effect = [
                np.array([1.0, 0.0, 0.0]),  # item
                np.array([0.95, 0.05, 0.0]),  # existing_1 (very similar)
            ]
            mock_explain.return_value = "Similar to existing cancer studies"
            
            score = score_novelty(item, existing_items)
            
            assert isinstance(score, NoveltyScore)
            assert score.novelty_score < 0.2  # Low novelty due to high similarity
            assert score.max_similarity > 0.8  # High similarity to existing
            assert score.most_similar_item_id == "existing_1"

    def test_batch_score_comprehensive(self):
        """Test batch scoring with both relevance and novelty."""
        items = [
            {
                "id": "item_1",
                "title": "Relevant novel study",
                "description": "New approach to cancer research"
            },
            {
                "id": "item_2",
                "title": "Less relevant study", 
                "description": "Standard analysis"
            }
        ]
        
        context = {"diseases": ["cancer"], "assay_types": ["RNA-seq"]}
        
        with patch('amprenta_rag.analysis.relevance_scoring.score_relevance') as mock_relevance, \
             patch('amprenta_rag.analysis.relevance_scoring.score_novelty') as mock_novelty:
            
            # Mock relevance scores
            mock_relevance.side_effect = [
                RelevanceScore("item_1", 0.9, 0.9, 0.8, 0.9, "High relevance", 1.0),
                RelevanceScore("item_2", 0.5, 0.6, 0.4, 0.5, "Medium relevance", 1.0)
            ]
            
            # Mock novelty scores  
            mock_novelty.side_effect = [
                NoveltyScore("item_1", 0.8, 0.2, "item_2", "Novel approach", 1.0),
                NoveltyScore("item_2", 0.3, 0.7, "item_1", "Similar to existing", 1.0)
            ]
            
            # Use correct parameter names: score_relevance_flag, score_novelty_flag
            scored_items = batch_score(items, context, score_relevance_flag=True, score_novelty_flag=True)
            
            assert len(scored_items) == 2
            assert all(isinstance(item, ScoredItem) for item in scored_items)
            
            # Check first item (should be ranked higher)
            item1 = next(item for item in scored_items if item.item_id == "item_1")
            assert item1.relevance_score.overall_score == 0.9
            assert item1.novelty_score.novelty_score == 0.8

    def test_batch_score_relevance_only(self):
        """Test batch scoring with relevance only."""
        items = [{"id": "item_1", "title": "Test study"}]
        context = {"diseases": ["cancer"]}
        
        with patch('amprenta_rag.analysis.relevance_scoring.score_relevance') as mock_relevance:
            mock_relevance.return_value = RelevanceScore("item_1", 0.7, 0.7, 0.6, 0.8, "Good relevance", 1.0)
            
            # Use correct parameter names: score_relevance_flag, score_novelty_flag
            scored_items = batch_score(items, context, score_relevance_flag=True, score_novelty_flag=False)
            
            assert len(scored_items) == 1
            assert scored_items[0].relevance_score.overall_score == 0.7
            assert scored_items[0].novelty_score is None

    def test_caching_functionality(self):
        """Test score caching and retrieval."""
        cache_key = "test_relevance_item123_context456"
        score_data = {"overall_score": 0.8, "explanation": "Test score"}
        
        # Test caching
        _set_cache(cache_key, score_data)
        
        # Test retrieval
        cached_score = _get_from_cache(cache_key)
        assert cached_score is not None
        assert cached_score["overall_score"] == 0.8


class TestDataClasses:
    """Test data class functionality."""

    def test_relevance_score_creation(self):
        """Test RelevanceScore data class creation."""
        score = RelevanceScore(
            item_id="test_123",
            overall_score=0.85,
            disease_match=0.9,
            target_overlap=0.8,
            data_quality=0.85,
            explanation="High relevance score",
            processing_time_seconds=2.5,
            cached=False
        )
        
        assert score.item_id == "test_123"
        assert score.overall_score == 0.85
        assert score.cached is False

    def test_novelty_score_creation(self):
        """Test NoveltyScore data class creation."""
        score = NoveltyScore(
            item_id="test_123",
            novelty_score=0.75,
            max_similarity=0.25,
            most_similar_item_id="similar_456",
            explanation="Moderately novel",
            processing_time_seconds=1.8,
            cached=True
        )
        
        assert score.item_id == "test_123"
        assert score.novelty_score == 0.75
        assert score.most_similar_item_id == "similar_456"
        assert score.cached is True

    def test_scored_item_creation(self):
        """Test ScoredItem data class creation."""
        relevance_score = RelevanceScore("test_123", 0.8, 0.8, 0.7, 0.9, "Relevant", 1.0)
        novelty_score = NoveltyScore("test_123", 0.6, 0.4, "other_item", "Somewhat novel", 1.0)
        
        # ScoredItem only has item_id, relevance_score, novelty_score
        scored_item = ScoredItem(
            item_id="test_123",
            relevance_score=relevance_score,
            novelty_score=novelty_score,
        )
        
        assert scored_item.item_id == "test_123"
        assert scored_item.relevance_score.overall_score == 0.8
        assert scored_item.novelty_score.novelty_score == 0.6


@pytest.mark.integration
class TestRelevanceScoringIntegration:
    """Integration tests for relevance scoring (requires OpenAI API)."""

    @pytest.mark.skip(reason="Requires OpenAI API access")
    def test_real_relevance_scoring(self):
        """Test real relevance scoring with OpenAI API."""
        item = {
            "id": "test_dataset",
            "title": "Breast cancer gene expression analysis",
            "description": "RNA-seq study of breast cancer patients",
            "species": "human"
        }
        
        context = {
            "diseases": ["breast cancer"],
            "assay_types": ["RNA-seq"],
            "species": ["human"]
        }
        
        score = score_relevance(item, context)
        
        assert isinstance(score, RelevanceScore)
        assert 0.0 <= score.overall_score <= 1.0
        assert len(score.explanation) > 0

    @pytest.mark.skip(reason="Requires OpenAI API access")
    def test_real_novelty_scoring(self):
        """Test real novelty scoring with OpenAI API."""
        item = {
            "id": "new_item",
            "title": "Novel cancer biomarker discovery",
            "description": "Identification of new prognostic markers"
        }
        
        existing_items = [
            {
                "id": "existing_1",
                "title": "Traditional cancer markers",
                "description": "Analysis of known cancer biomarkers"
            }
        ]
        
        score = score_novelty(item, existing_items)
        
        assert isinstance(score, NoveltyScore)
        assert 0.0 <= score.novelty_score <= 1.0
        assert 0.0 <= score.max_similarity <= 1.0
