"""
Unit tests for scoring API endpoints.

Tests relevance and novelty scoring endpoints.
Note: These endpoints use lazy imports, so we test with real implementations.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

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
