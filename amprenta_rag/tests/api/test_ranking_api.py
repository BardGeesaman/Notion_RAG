"""Tests for compound ranking API endpoints."""

from __future__ import annotations

import asyncio
from types import SimpleNamespace
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app


pytest.importorskip("fastapi")

client = TestClient(app)


def test_get_presets() -> None:
    """Test GET /api/ranking/presets endpoint."""
    resp = client.get("/api/ranking/presets")
    assert resp.status_code == 200
    
    data = resp.json()
    assert isinstance(data, list)
    assert len(data) >= 4  # At least the 4 built-in presets
    
    # Check structure of first preset
    preset = data[0]
    assert "name" in preset
    assert "weights" in preset
    assert "description" in preset
    assert isinstance(preset["weights"], dict)
    
    # Check that balanced preset exists
    preset_names = [p["name"] for p in data]
    assert "balanced" in preset_names


def test_score_endpoint_with_preset(monkeypatch) -> None:
    """Test POST /api/ranking/score endpoint with preset weights."""
    from amprenta_rag.api.routers import ranking as ranking_router
    
    # Mock rank_compounds service
    mock_rankings = [
        SimpleNamespace(
            compound_id=str(uuid4()),
            smiles="CCO",
            objectives=[
                SimpleNamespace(name="potency", raw_value=100.0, normalized=0.8, weight=0.4, confidence=None),
                SimpleNamespace(name="herg", raw_value=0.2, normalized=0.8, weight=0.3, confidence=0.1),
            ],
            weighted_score=0.8,
            pareto_rank=1,
            rank=1
        )
    ]
    
    def mock_rank_compounds(compound_ids, weights, db, include_pareto=True):
        return mock_rankings
    
    monkeypatch.setattr(ranking_router, "rank_compounds", mock_rank_compounds)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[ranking_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/ranking/score",
            json={
                "compound_ids": [str(uuid4())],
                "preset": "balanced",
                "include_pareto": True
            }
        )
        assert resp.status_code == 200
        
        data = resp.json()
        assert "rankings" in data
        assert "pareto_front" in data
        assert "total_compounds" in data
        assert "skipped_compounds" in data
        
        assert len(data["rankings"]) == 1
        assert len(data["pareto_front"]) == 1  # Mock ranking has pareto_rank=1
        assert data["total_compounds"] == 1
        
        ranking = data["rankings"][0]
        assert "compound_id" in ranking
        assert "smiles" in ranking
        assert "objectives" in ranking
        assert "weighted_score" in ranking
        assert "pareto_rank" in ranking
        assert "rank" in ranking
        
    finally:
        del app.dependency_overrides[ranking_router.get_db]


def test_score_endpoint_with_custom_weights(monkeypatch) -> None:
    """Test POST /api/ranking/score endpoint with custom weights."""
    from amprenta_rag.api.routers import ranking as ranking_router
    
    # Mock rank_compounds service
    def mock_rank_compounds(compound_ids, weights, db, include_pareto=True):
        # Verify custom weights were passed
        assert weights["potency"] == 0.6
        assert weights["herg"] == 0.4
        return []
    
    monkeypatch.setattr(ranking_router, "rank_compounds", mock_rank_compounds)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[ranking_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/ranking/score",
            json={
                "compound_ids": [str(uuid4())],
                "weights": {"potency": 0.6, "herg": 0.4},
                "include_pareto": False
            }
        )
        assert resp.status_code == 200
        
    finally:
        del app.dependency_overrides[ranking_router.get_db]


def test_pareto_endpoint(monkeypatch) -> None:
    """Test POST /api/ranking/pareto endpoint."""
    from amprenta_rag.api.routers import ranking as ranking_router
    
    # Mock rank_compounds service
    comp1_id = str(uuid4())
    comp2_id = str(uuid4())
    mock_rankings = [
        SimpleNamespace(
            compound_id=comp1_id,
            smiles="CCO",
            objectives=[
                SimpleNamespace(name="potency", normalized=0.8),
                SimpleNamespace(name="herg", normalized=0.7),
                SimpleNamespace(name="alerts", normalized=0.9),
            ],
            pareto_rank=1
        ),
        SimpleNamespace(
            compound_id=comp2_id, 
            smiles="CCC",
            objectives=[
                SimpleNamespace(name="potency", normalized=0.6),
                SimpleNamespace(name="herg", normalized=0.5),
                SimpleNamespace(name="alerts", normalized=0.8),
            ],
            pareto_rank=2
        )
    ]
    
    def mock_rank_compounds(compound_ids, weights, db, include_pareto=True):
        return mock_rankings
    
    monkeypatch.setattr(ranking_router, "rank_compounds", mock_rank_compounds)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[ranking_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/ranking/pareto",
            json={
                "compound_ids": [comp1_id, comp2_id],
                "x_objective": "potency",
                "y_objective": "liability_aggregate"
            }
        )
        assert resp.status_code == 200
        
        data = resp.json()
        assert "points" in data
        assert "x_label" in data
        assert "y_label" in data
        assert "frontier_ids" in data
        
        assert len(data["points"]) == 2
        assert comp1_id in data["frontier_ids"]  # pareto_rank=1
        assert comp2_id not in data["frontier_ids"]  # pareto_rank=2
        
        # Check point structure
        point = data["points"][0]
        assert "compound_id" in point
        assert "smiles" in point
        assert "x_value" in point
        assert "y_value" in point
        assert "pareto_rank" in point
        assert "is_frontier" in point
        
    finally:
        del app.dependency_overrides[ranking_router.get_db]


def test_score_endpoint_empty_compounds() -> None:
    """Test score endpoint with empty compound list."""
    resp = client.post(
        "/api/ranking/score",
        json={"compound_ids": [], "preset": "balanced"}
    )
    assert resp.status_code == 400
    assert "No compound IDs provided" in resp.json()["detail"]


def test_score_endpoint_exceeds_limit() -> None:
    """Test score endpoint with too many compounds."""
    # Generate more than MAX_COMPOUNDS (500) UUIDs
    many_ids = [str(uuid4()) for _ in range(501)]
    
    resp = client.post(
        "/api/ranking/score",
        json={"compound_ids": many_ids, "preset": "balanced"}
    )
    assert resp.status_code == 400
    assert "Too many compounds" in resp.json()["detail"]
    assert "Maximum allowed: 500" in resp.json()["detail"]


def test_score_endpoint_invalid_uuid() -> None:
    """Test score endpoint with invalid UUID format."""
    resp = client.post(
        "/api/ranking/score",
        json={"compound_ids": ["not-a-uuid"], "preset": "balanced"}
    )
    assert resp.status_code == 400
    assert "Invalid UUID format" in resp.json()["detail"]


def test_pareto_endpoint_invalid_uuid() -> None:
    """Test pareto endpoint with invalid UUID format."""
    resp = client.post(
        "/api/ranking/pareto",
        json={
            "compound_ids": ["invalid-uuid"],
            "x_objective": "potency",
            "y_objective": "herg"
        }
    )
    assert resp.status_code == 400
    assert "Invalid UUID format" in resp.json()["detail"]


def test_score_endpoint_unknown_preset(monkeypatch) -> None:
    """Test score endpoint with unknown preset defaults to balanced."""
    from amprenta_rag.api.routers import ranking as ranking_router
    
    # Mock rank_compounds service to check weights
    def mock_rank_compounds(compound_ids, weights, db, include_pareto=True):
        # Should default to balanced preset weights
        assert "potency" in weights
        assert "herg" in weights
        return []
    
    monkeypatch.setattr(ranking_router, "rank_compounds", mock_rank_compounds)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[ranking_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/ranking/score",
            json={
                "compound_ids": [str(uuid4())],
                "preset": "unknown_preset"
            }
        )
        assert resp.status_code == 200
        
    finally:
        del app.dependency_overrides[ranking_router.get_db]


def test_pareto_endpoint_specific_objective(monkeypatch) -> None:
    """Test pareto endpoint with specific Y objective (not aggregate)."""
    from amprenta_rag.api.routers import ranking as ranking_router
    
    # Mock rank_compounds service
    comp1_id = str(uuid4())
    mock_rankings = [
        SimpleNamespace(
            compound_id=comp1_id,
            smiles="CCO",
            objectives=[
                SimpleNamespace(name="potency", normalized=0.8),
                SimpleNamespace(name="herg", normalized=0.7),
            ],
            pareto_rank=1
        )
    ]
    
    def mock_rank_compounds(compound_ids, weights, db, include_pareto=True):
        return mock_rankings
    
    monkeypatch.setattr(ranking_router, "rank_compounds", mock_rank_compounds)
    
    # Mock database dependency
    def mock_get_db():
        yield SimpleNamespace()
    
    app.dependency_overrides[ranking_router.get_db] = mock_get_db
    
    try:
        resp = client.post(
            "/api/ranking/pareto",
            json={
                "compound_ids": [comp1_id],
                "x_objective": "potency",
                "y_objective": "herg"  # Specific objective, not aggregate
            }
        )
        assert resp.status_code == 200
        
        data = resp.json()
        assert data["y_label"] == "Herg"  # Should format objective name
        
    finally:
        del app.dependency_overrides[ranking_router.get_db]


class TestAsyncLLMRanking:
    """Test async execution of LLM-based ranking endpoints."""

    @pytest.mark.asyncio
    async def test_rank_async(self):
        """Test rank endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.ranking._sync_rank_items') as mock_rank:
            from amprenta_rag.api.schemas import RankingItem
            
            # Mock the ranking function
            mock_result = MagicMock()
            mock_result.ranked_items = [
                RankingItem(
                    item_id="item_1",
                    rank=1,
                    score=0.95,
                    explanation="Highest relevance"
                ),
                RankingItem(
                    item_id="item_2",
                    rank=2,
                    score=0.82,
                    explanation="Good relevance"
                )
            ]
            mock_result.criteria_used = "Relevance to ALS research"
            mock_result.explanation = "Ranked by disease relevance and data quality"
            mock_result.processing_time_seconds = 0.5
            mock_result.cached = False
            mock_rank.return_value = mock_result
            
            from amprenta_rag.api.routers.ranking import rank_items_endpoint
            from amprenta_rag.api.schemas import RankRequest, ScoringItemRequest, ScoringContextRequest
            
            request = RankRequest(
                items=[
                    ScoringItemRequest(
                        id="item_1",
                        title="ALS RNA-seq Dataset",
                        description="Comprehensive RNA-seq data from ALS patients",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=50,
                    ),
                    ScoringItemRequest(
                        id="item_2",
                        title="ALS Proteomics Dataset",
                        description="Proteomics analysis of ALS samples",
                        species="human",
                        assay_type="proteomics",
                        sample_count=30,
                    ),
                ],
                criteria="Relevance to ALS research",
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=["SOD1"],
                    species=["human"],
                    assay_types=["RNA-seq", "proteomics"],
                    min_sample_size=20,
                ),
            )
            
            result = await rank_items_endpoint(request)
            
            # Verify async execution and result
            assert len(result.ranked_items) == 2
            assert result.ranked_items[0].item_id == "item_1"
            assert result.ranked_items[0].rank == 1
            assert result.ranked_items[1].item_id == "item_2"
            assert result.ranked_items[1].rank == 2
            assert result.criteria_used == "Relevance to ALS research"
            mock_rank.assert_called_once()

    @pytest.mark.asyncio
    async def test_rerank_async(self):
        """Test rerank endpoint executes asynchronously."""
        with patch('amprenta_rag.api.routers.ranking._sync_rerank_items') as mock_rerank:
            from amprenta_rag.api.schemas import RankingItem
            
            # Mock the re-ranking function
            mock_result = MagicMock()
            mock_result.ranked_items = [
                RankingItem(
                    item_id="item_2",
                    rank=1,
                    score=0.88,
                    explanation="Re-ranked higher due to sample size"
                ),
                RankingItem(
                    item_id="item_1",
                    rank=2,
                    score=0.85,
                    explanation="Re-ranked lower due to new criteria"
                )
            ]
            mock_result.criteria_used = "Sample size priority"
            mock_result.explanation = "Re-ranked prioritizing larger sample sizes"
            mock_result.processing_time_seconds = 0.3
            mock_result.cached = False
            mock_rerank.return_value = mock_result
            
            from amprenta_rag.api.routers.ranking import rerank_items_endpoint
            from amprenta_rag.api.schemas import RerankRequest, ScoringItemRequest, ScoringContextRequest
            
            request = RerankRequest(
                items=[
                    ScoringItemRequest(
                        id="item_1",
                        title="Small Dataset",
                        description="Dataset with fewer samples",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=20,
                    ),
                    ScoringItemRequest(
                        id="item_2",
                        title="Large Dataset",
                        description="Dataset with more samples",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=100,
                    ),
                ],
                previous_ranking=["item_1", "item_2"],
                criteria="Sample size priority",
                context=ScoringContextRequest(
                    diseases=["ALS"],
                    targets=[],
                    species=["human"],
                    assay_types=["RNA-seq"],
                    min_sample_size=50,
                ),
            )
            
            result = await rerank_items_endpoint(request)
            
            # Verify async execution and result
            assert len(result.ranked_items) == 2
            assert result.ranked_items[0].item_id == "item_2"  # Re-ranked to top
            assert result.ranked_items[0].rank == 1
            assert result.ranked_items[1].item_id == "item_1"  # Re-ranked to bottom
            assert result.ranked_items[1].rank == 2
            assert result.criteria_used == "Sample size priority"
            mock_rerank.assert_called_once()

    @pytest.mark.asyncio
    async def test_ranking_concurrent(self):
        """Test multiple simultaneous ranking requests."""
        with patch('amprenta_rag.api.routers.ranking._sync_rank_items') as mock_rank:
            from amprenta_rag.api.schemas import RankingItem
            
            # Mock function to return different results for different calls
            def mock_rank_func(items, criteria, context):
                mock_result = MagicMock()
                mock_result.ranked_items = [
                    RankingItem(
                        item_id=items[0]["id"],
                        rank=1,
                        score=0.9,
                        explanation=f"Top ranked for {criteria}"
                    )
                ]
                mock_result.criteria_used = criteria
                mock_result.explanation = f"Ranked using {criteria}"
                mock_result.processing_time_seconds = 0.2
                mock_result.cached = False
                return mock_result
            
            mock_rank.side_effect = mock_rank_func
            
            from amprenta_rag.api.routers.ranking import rank_items_endpoint
            from amprenta_rag.api.schemas import RankRequest, ScoringItemRequest, ScoringContextRequest
            
            async def make_ranking_request(item_id: str, criteria: str):
                request = RankRequest(
                    items=[
                        ScoringItemRequest(
                            id=item_id,
                            title=f"Dataset {item_id}",
                            description=f"Test dataset {item_id}",
                            species="human",
                            assay_type="RNA-seq",
                            sample_count=25,
                        )
                    ],
                    criteria=criteria,
                    context=ScoringContextRequest(
                        diseases=["ALS"],
                        targets=[],
                        species=["human"],
                        assay_types=["RNA-seq"],
                        min_sample_size=10,
                    ),
                )
                return await rank_items_endpoint(request)
            
            # Make 3 concurrent requests
            tasks = [
                make_ranking_request("concurrent_1", "criteria_1"),
                make_ranking_request("concurrent_2", "criteria_2"),
                make_ranking_request("concurrent_3", "criteria_3"),
            ]
            
            results = await asyncio.gather(*tasks)
            
            # All requests should succeed
            assert len(results) == 3
            for i, result in enumerate(results, 1):
                assert result.ranked_items[0].item_id == f"concurrent_{i}"
                assert result.criteria_used == f"criteria_{i}"
                assert f"criteria_{i}" in result.explanation
            
            # Verify all calls were made
            assert mock_rank.call_count == 3

    @pytest.mark.asyncio
    async def test_ranking_error_handling(self):
        """Test exception propagation in async ranking."""
        with patch('amprenta_rag.api.routers.ranking._sync_rank_items') as mock_rank:
            # Mock function to raise an exception
            mock_rank.side_effect = Exception("LLM service unavailable")
            
            from amprenta_rag.api.routers.ranking import rank_items_endpoint
            from amprenta_rag.api.schemas import RankRequest, ScoringItemRequest
            
            request = RankRequest(
                items=[
                    ScoringItemRequest(
                        id="error_test",
                        title="Error Test Dataset",
                        description="Dataset for testing error handling",
                        species="human",
                        assay_type="RNA-seq",
                        sample_count=20,
                    )
                ],
                criteria="Test criteria",
                context=None,
            )
            
            # Should raise HTTPException
            with pytest.raises(Exception) as exc_info:
                await rank_items_endpoint(request)
            
            # Verify the exception is properly handled
            # (In a real FastAPI app, this would be converted to HTTPException)
            assert "LLM service unavailable" in str(exc_info.value)
            mock_rank.assert_called_once()
