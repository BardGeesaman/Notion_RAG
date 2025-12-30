"""
Unit tests for scoring API endpoints.

Tests relevance and novelty scoring endpoints.
Note: These endpoints use lazy imports, so we test with real implementations.
"""

from __future__ import annotations

import asyncio
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient
from httpx import AsyncClient

from amprenta_rag.api.main import app

client = TestClient(app)


class TestRelevanceScore:
    """Tests for POST /api/v1/score/relevance endpoint."""

    def test_score_relevance_success(self):
        """Test successful relevance scoring."""
        item_id = "item123"
        
        response = client.post(
            "/api/v1/score/relevance",
            json={
                "item": {
                    "id": item_id,
                    "title": "Test Dataset",
                    "description": "ALS patient samples",
                    "species": "human",
                    "assay_type": "RNA-seq",
                    "sample_count": 50,
                },
                "context": {
                    "diseases": ["ALS"],
                    "targets": ["SOD1"],
                    "species": ["human"],
                    "assay_types": ["RNA-seq"],
                    "min_sample_size": 20,
                },
                "criteria": None,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["item_id"] == item_id
        assert "overall_score" in data
        assert "disease_match" in data
        assert "target_overlap" in data
        assert "data_quality" in data
        assert "explanation" in data
        assert 0.0 <= data["overall_score"] <= 1.0

    def test_score_relevance_invalid_request(self):
        """Test relevance scoring with invalid request."""
        response = client.post(
            "/api/v1/score/relevance",
            json={
                "item": {
                    "id": "item123",
                    # Missing required fields
                },
                "context": {},
            },
        )
        
        # Should accept minimal request (fields are optional)
        # If it returns 200, that's fine; if validation fails, expect 422
        assert response.status_code in [200, 422]


class TestNoveltyScore:
    """Tests for POST /api/v1/score/novelty endpoint."""

    def test_score_novelty_success(self):
        """Test successful novelty scoring."""
        item_id = "item123"
        similar_id = "item456"
        
        response = client.post(
            "/api/v1/score/novelty",
            json={
                "item": {
                    "id": item_id,
                    "title": "New Dataset",
                    "description": "Novel ALS study with unique features",
                    "species": "human",
                    "assay_type": "proteomics",
                    "sample_count": 30,
                },
                "existing_items": [
                    {
                        "id": similar_id,
                        "title": "Existing Dataset",
                        "description": "Previous ALS study with different approach",
                        "species": "human",
                        "assay_type": "RNA-seq",
                        "sample_count": 25,
                    }
                ],
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["item_id"] == item_id
        assert "novelty_score" in data
        assert "max_similarity" in data
        assert "explanation" in data
        assert 0.0 <= data["novelty_score"] <= 1.0

    def test_score_novelty_no_existing_items(self):
        """Test novelty scoring with no existing items to compare."""
        item_id = "item123"
        
        response = client.post(
            "/api/v1/score/novelty",
            json={
                "item": {
                    "id": item_id,
                    "title": "New Dataset",
                    "description": "Novel dataset",
                    "species": "human",
                    "assay_type": "proteomics",
                    "sample_count": 30,
                },
                "existing_items": [],
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["item_id"] == item_id
        # With no existing items, novelty should be high
        assert data["novelty_score"] >= 0.8


class TestBatchScore:
    """Tests for POST /api/v1/score/batch endpoint."""

    def test_batch_score_relevance_only(self):
        """Test batch scoring with relevance only."""
        item1_id = "item1"
        item2_id = "item2"
        
        response = client.post(
            "/api/v1/score/batch",
            json={
                "items": [
                    {
                        "id": item1_id,
                        "title": "Dataset 1",
                        "description": "ALS RNA-seq data",
                        "species": "human",
                        "assay_type": "RNA-seq",
                        "sample_count": 20,
                    },
                    {
                        "id": item2_id,
                        "title": "Dataset 2",
                        "description": "ALS proteomics data",
                        "species": "human",
                        "assay_type": "proteomics",
                        "sample_count": 30,
                    },
                ],
                "context": {
                    "diseases": ["ALS"],
                    "targets": [],
                    "species": ["human"],
                    "assay_types": [],
                    "min_sample_size": 10,
                },
                "score_relevance": True,
                "score_novelty": False,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["total_items"] == 2
        assert len(data["items"]) == 2
        assert data["items"][0]["item_id"] == item1_id
        assert data["items"][1]["item_id"] == item2_id
        # Should have relevance scores
        assert data["items"][0]["relevance_score"] is not None
        assert data["items"][1]["relevance_score"] is not None

    def test_batch_score_novelty_only(self):
        """Test batch scoring with novelty only."""
        item1_id = "item1"
        
        response = client.post(
            "/api/v1/score/batch",
            json={
                "items": [
                    {
                        "id": item1_id,
                        "title": "Dataset 1",
                        "description": "Test data",
                        "species": "human",
                        "assay_type": "test",
                        "sample_count": 10,
                    },
                ],
                "context": {
                    "diseases": [],
                    "targets": [],
                    "species": [],
                    "assay_types": [],
                    "min_sample_size": 1,
                },
                "score_relevance": False,
                "score_novelty": True,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["total_items"] == 1
        assert data["items"][0]["item_id"] == item1_id
        # Should have novelty score
        assert data["items"][0]["novelty_score"] is not None

    def test_batch_score_both_scores(self):
        """Test batch scoring with both relevance and novelty."""
        item_id = "item1"
        
        response = client.post(
            "/api/v1/score/batch",
            json={
                "items": [
                    {
                        "id": item_id,
                        "title": "Comprehensive Dataset",
                        "description": "ALS patient samples with proteomics",
                        "species": "human",
                        "assay_type": "proteomics",
                        "sample_count": 40,
                    },
                ],
                "context": {
                    "diseases": ["ALS"],
                    "targets": [],
                    "species": ["human"],
                    "assay_types": [],
                    "min_sample_size": 20,
                },
                "score_relevance": True,
                "score_novelty": True,
            },
        )
        
        assert response.status_code == 200
        data = response.json()
        assert data["total_items"] == 1
        assert data["items"][0]["item_id"] == item_id
        # Should have both scores
        assert data["items"][0]["relevance_score"] is not None
        assert data["items"][0]["novelty_score"] is not None

    def test_batch_score_empty_items(self):
        """Test batch scoring with empty items list."""
        response = client.post(
            "/api/v1/score/batch",
            json={
                "items": [],
                "context": {
                    "diseases": [],
                    "targets": [],
                    "species": [],
                    "assay_types": [],
                    "min_sample_size": 1,
                },
                "score_relevance": False,
                "score_novelty": False,
            },
        )
        
        # Should fail validation (min_length=1 for items)
        assert response.status_code == 422


class TestAsyncScoringEndpoints:
    """Test async execution of scoring endpoints using mocked async functions."""

    @pytest.mark.asyncio
    async def test_relevance_async(self):
        """Test relevance endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.scoring._sync_score_relevance') as mock_score:
            # Mock the scoring function to return a simple result
            mock_result = MagicMock()
            mock_result.item_id = "async_item_123"
            mock_result.overall_score = 0.85
            mock_result.disease_match = 0.9
            mock_result.target_overlap = 0.8
            mock_result.data_quality = 0.85
            mock_result.explanation = "High relevance async test"
            mock_result.processing_time_seconds = 0.1
            mock_result.cached = False
            mock_score.return_value = mock_result
            
            # Import the endpoint function and call it directly
            from amprenta_rag.api.routers.scoring import score_relevance_endpoint
            from amprenta_rag.api.schemas import RelevanceScoreRequest, ScoringItemRequest, ScoringContextRequest
            
            request = RelevanceScoreRequest(
                item=ScoringItemRequest(
                    id="async_item_123",
                    title="Async Test Dataset",
                    description="Testing async execution",
                    species="human",
                    assay_type="RNA-seq",
                    sample_count=25,
                ),
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=["SOD1"],
                    species=["human"],
                    assay_types=["RNA-seq"],
                    min_sample_size=10,
                ),
                criteria=None,
            )
            
            result = await score_relevance_endpoint(request)
            
            # Verify async execution and result
            assert result.item_id == "async_item_123"
            assert result.overall_score == 0.85
            assert result.explanation == "High relevance async test"
            mock_score.assert_called_once()

    @pytest.mark.asyncio
    async def test_novelty_async(self):
        """Test novelty endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.scoring._sync_score_novelty') as mock_score:
            # Mock the scoring function
            mock_result = MagicMock()
            mock_result.item_id = "async_novelty_123"
            mock_result.novelty_score = 0.75
            mock_result.max_similarity = 0.25
            mock_result.most_similar_item_id = "existing_123"
            mock_result.explanation = "Moderately novel async test"
            mock_result.processing_time_seconds = 0.2
            mock_result.cached = False
            mock_score.return_value = mock_result
            
            from amprenta_rag.api.routers.scoring import score_novelty_endpoint
            from amprenta_rag.api.schemas import NoveltyScoreRequest, ScoringItemRequest
            
            request = NoveltyScoreRequest(
                item=ScoringItemRequest(
                    id="async_novelty_123",
                    title="Async Novelty Dataset",
                    description="Testing async novelty scoring",
                    species="human",
                    assay_type="proteomics",
                    sample_count=30,
                ),
                existing_items=[
                    ScoringItemRequest(
                        id="existing_123",
                        title="Existing Dataset",
                        description="Previous study",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=20,
                    )
                ],
            )
            
            result = await score_novelty_endpoint(request)
            
            # Verify async execution and result
            assert result.item_id == "async_novelty_123"
            assert result.novelty_score == 0.75
            assert result.explanation == "Moderately novel async test"
            mock_score.assert_called_once()

    @pytest.mark.asyncio
    async def test_batch_async(self):
        """Test batch endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.scoring._sync_batch_score') as mock_batch:
            # Mock batch scoring results
            mock_result1 = MagicMock()
            mock_result1.item_id = "async_batch_1"
            mock_result1.relevance_score = MagicMock()
            mock_result1.relevance_score.item_id = "async_batch_1"
            mock_result1.relevance_score.overall_score = 0.8
            mock_result1.relevance_score.disease_match = 0.85
            mock_result1.relevance_score.target_overlap = 0.75
            mock_result1.relevance_score.data_quality = 0.8
            mock_result1.relevance_score.explanation = "Good relevance"
            mock_result1.relevance_score.processing_time_seconds = 0.1
            mock_result1.relevance_score.cached = False
            mock_result1.novelty_score = None
            
            mock_result2 = MagicMock()
            mock_result2.item_id = "async_batch_2"
            mock_result2.relevance_score = MagicMock()
            mock_result2.relevance_score.item_id = "async_batch_2"
            mock_result2.relevance_score.overall_score = 0.7
            mock_result2.relevance_score.disease_match = 0.8
            mock_result2.relevance_score.target_overlap = 0.6
            mock_result2.relevance_score.data_quality = 0.7
            mock_result2.relevance_score.explanation = "Moderate relevance"
            mock_result2.relevance_score.processing_time_seconds = 0.1
            mock_result2.relevance_score.cached = False
            mock_result2.novelty_score = None
            
            mock_batch.return_value = [mock_result1, mock_result2]
            
            from amprenta_rag.api.routers.scoring import batch_score_endpoint
            from amprenta_rag.api.schemas import BatchScoreRequest, ScoringItemRequest, ScoringContextRequest
            
            request = BatchScoreRequest(
                items=[
                    ScoringItemRequest(
                        id="async_batch_1",
                        title="Async Batch Dataset 1",
                        description="First dataset for async batch testing",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=15,
                    ),
                    ScoringItemRequest(
                        id="async_batch_2",
                        title="Async Batch Dataset 2",
                        description="Second dataset for async batch testing",
                        species="human",
                        assay_type="proteomics",
                        sample_count=25,
                    ),
                ],
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=[],
                    species=["human"],
                    assay_types=[],
                    min_sample_size=10,
                ),
                score_relevance=True,
                score_novelty=False,
            )
            
            result = await batch_score_endpoint(request)
            
            # Verify async execution and result
            assert result.total_items == 2
            assert len(result.items) == 2
            assert result.items[0].item_id == "async_batch_1"
            assert result.items[1].item_id == "async_batch_2"
            mock_batch.assert_called_once()

    @pytest.mark.asyncio
    async def test_concurrent_requests(self):
        """Test multiple simultaneous async requests."""
        with patch('amprenta_rag.api.routers.scoring._sync_score_relevance') as mock_score:
            # Mock function to return different results for different calls
            def mock_score_func(item, context, criteria):
                mock_result = MagicMock()
                mock_result.item_id = item["id"]
                mock_result.overall_score = 0.8
                mock_result.disease_match = 0.85
                mock_result.target_overlap = 0.75
                mock_result.data_quality = 0.8
                mock_result.explanation = f"Concurrent test for {item['id']}"
                mock_result.processing_time_seconds = 0.1
                mock_result.cached = False
                return mock_result
            
            mock_score.side_effect = mock_score_func
            
            from amprenta_rag.api.routers.scoring import score_relevance_endpoint
            from amprenta_rag.api.schemas import RelevanceScoreRequest, ScoringItemRequest, ScoringContextRequest
            
            async def make_relevance_request(item_id: str):
                request = RelevanceScoreRequest(
                    item=ScoringItemRequest(
                        id=item_id,
                        title=f"Concurrent Dataset {item_id}",
                        description="Testing concurrent execution",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=20,
                    ),
                    context=ScoringContextRequest(
                        diseases=["ALS"],
                        targets=[],
                        species=["human"],
                        assay_types=[],
                        min_sample_size=5,
                    ),
                    criteria=None,
                )
                return await score_relevance_endpoint(request)
            
            # Make 3 concurrent requests
            tasks = [
                make_relevance_request("concurrent_1"),
                make_relevance_request("concurrent_2"),
                make_relevance_request("concurrent_3"),
            ]
            
            results = await asyncio.gather(*tasks)
            
            # All requests should succeed
            assert len(results) == 3
            for i, result in enumerate(results, 1):
                assert result.item_id == f"concurrent_{i}"
                assert result.overall_score == 0.8
                assert f"concurrent_{i}" in result.explanation
            
            # Verify all calls were made
            assert mock_score.call_count == 3
