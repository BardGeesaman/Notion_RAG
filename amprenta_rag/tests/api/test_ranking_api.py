"""Tests for compound ranking API endpoints."""

from __future__ import annotations

from types import SimpleNamespace
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
