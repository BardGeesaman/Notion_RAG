"""Unit tests for Bayesian optimization module."""

from __future__ import annotations

import pytest


def test_recommend_next_compounds_synthetic():
    """Test BO with synthetic compound data."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")

    from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds

    # 5 tested compounds with 3D features
    tested = [
        {"features": [0.1, 0.2, 0.3], "activity": 0.5},
        {"features": [0.2, 0.3, 0.4], "activity": 0.6},
        {"features": [0.3, 0.4, 0.5], "activity": 0.7},
        {"features": [0.4, 0.5, 0.6], "activity": 0.8},
        {"features": [0.5, 0.6, 0.7], "activity": 0.9},
    ]

    candidates = [
        {"id": "cand_1", "features": [0.6, 0.7, 0.8]},
        {"id": "cand_2", "features": [0.7, 0.8, 0.9]},
        {"id": "cand_3", "features": [0.15, 0.25, 0.35]},
    ]

    result = recommend_next_compounds(tested, candidates, batch_size=2)

    assert len(result) == 2
    assert all("id" in r for r in result)
    assert all("acquisition" in r for r in result)
    assert all("pred_mean" in r for r in result)
    assert all("pred_std" in r for r in result)
    # Should be sorted by acquisition (descending)
    assert result[0]["acquisition"] >= result[1]["acquisition"]


def test_recommend_next_compounds_empty_pool():
    """Empty candidate pool returns empty list."""
    from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds

    tested = [{"features": [0.1], "activity": 0.5} for _ in range(5)]
    result = recommend_next_compounds(tested, [], batch_size=5)
    assert result == []


def test_recommend_next_compounds_validation():
    """Test input validation."""
    from amprenta_rag.analysis.bayesian_optimization import recommend_next_compounds

    # Too few tested compounds
    with pytest.raises(ValueError, match="at least 5"):
        recommend_next_compounds([{"features": [0.1], "activity": 0.5}], [], batch_size=1)


def test_recommend_multi_objective_2_objectives():
    """Test multi-objective BO with 2 objectives."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
    
    tested = [
        {"features": [0.1, 0.2], "potency": 0.5, "herg": 0.8},
        {"features": [0.2, 0.3], "potency": 0.6, "herg": 0.7},
        {"features": [0.3, 0.4], "potency": 0.7, "herg": 0.6},
        {"features": [0.4, 0.5], "potency": 0.8, "herg": 0.5},
        {"features": [0.5, 0.6], "potency": 0.9, "herg": 0.4},
    ]
    
    candidates = [
        {"id": "c1", "features": [0.6, 0.7]},
        {"id": "c2", "features": [0.7, 0.8]},
    ]
    
    result = recommend_multi_objective(
        tested_compounds=tested,
        candidate_pool=candidates,
        objectives=["potency", "herg"],
        objective_directions=["maximize", "minimize"],
        batch_size=2,
    )
    
    assert "recommendations" in result
    assert "pareto_front" in result
    assert len(result["recommendations"]) == 2
    assert all("id" in r for r in result["recommendations"])
    assert all("acquisition" in r for r in result["recommendations"])
    assert all("pred_means" in r for r in result["recommendations"])
    assert all("pred_stds" in r for r in result["recommendations"])
    # Check predictions include both objectives
    for rec in result["recommendations"]:
        assert "potency" in rec["pred_means"]
        assert "herg" in rec["pred_means"]


def test_recommend_multi_objective_3_objectives():
    """Test multi-objective BO with 3 objectives."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
    
    tested = [
        {"features": [0.1, 0.2], "potency": 0.5, "herg": 0.8, "logs": -2.0},
        {"features": [0.2, 0.3], "potency": 0.6, "herg": 0.7, "logs": -1.5},
        {"features": [0.3, 0.4], "potency": 0.7, "herg": 0.6, "logs": -1.0},
        {"features": [0.4, 0.5], "potency": 0.8, "herg": 0.5, "logs": -0.8},
        {"features": [0.5, 0.6], "potency": 0.9, "herg": 0.4, "logs": -0.5},
    ]
    
    candidates = [
        {"id": "c1", "features": [0.6, 0.7]},
        {"id": "c2", "features": [0.7, 0.8]},
    ]
    
    result = recommend_multi_objective(
        tested_compounds=tested,
        candidate_pool=candidates,
        objectives=["potency", "herg", "logs"],
        objective_directions=["maximize", "minimize", "maximize"],
        batch_size=2,
    )
    
    assert len(result["recommendations"]) == 2
    # Check all 3 objectives in predictions
    for rec in result["recommendations"]:
        assert "potency" in rec["pred_means"]
        assert "herg" in rec["pred_means"]
        assert "logs" in rec["pred_means"]


def test_objective_directions_minimize():
    """Test that minimize direction works correctly."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
    
    tested = [
        {"features": [0.1], "obj1": 0.2, "obj2": 0.8},
        {"features": [0.2], "obj1": 0.3, "obj2": 0.7},
        {"features": [0.3], "obj1": 0.4, "obj2": 0.6},
        {"features": [0.4], "obj1": 0.5, "obj2": 0.5},
        {"features": [0.5], "obj1": 0.6, "obj2": 0.4},
    ]
    
    candidates = [{"id": "c1", "features": [0.6]}]
    
    # Should not raise error with minimize direction
    result = recommend_multi_objective(
        tested_compounds=tested,
        candidate_pool=candidates,
        objectives=["obj1", "obj2"],
        objective_directions=["maximize", "minimize"],
        batch_size=1,
    )
    
    assert len(result["recommendations"]) == 1


def test_validation_missing_objective():
    """Test validation fails when objective missing from tested data."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
    
    tested = [
        {"features": [0.1], "potency": 0.5},  # Missing herg
        {"features": [0.2], "potency": 0.6},
        {"features": [0.3], "potency": 0.7},
        {"features": [0.4], "potency": 0.8},
        {"features": [0.5], "potency": 0.9},
    ]
    
    candidates = [{"id": "c1", "features": [0.6]}]
    
    with pytest.raises(ValueError, match="missing from tested compound"):
        recommend_multi_objective(
            tested_compounds=tested,
            candidate_pool=candidates,
            objectives=["potency", "herg"],
            objective_directions=["maximize", "minimize"],
            batch_size=1,
        )


def test_validation_direction_mismatch():
    """Test validation fails when objectives and directions length mismatch."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import recommend_multi_objective
    
    tested = [
        {"features": [0.1], "potency": 0.5, "herg": 0.8},
        {"features": [0.2], "potency": 0.6, "herg": 0.7},
        {"features": [0.3], "potency": 0.7, "herg": 0.6},
        {"features": [0.4], "potency": 0.8, "herg": 0.5},
        {"features": [0.5], "potency": 0.9, "herg": 0.4},
    ]
    
    candidates = [{"id": "c1", "features": [0.6]}]
    
    with pytest.raises(ValueError, match="same length"):
        recommend_multi_objective(
            tested_compounds=tested,
            candidate_pool=candidates,
            objectives=["potency", "herg"],
            objective_directions=["maximize"],  # Wrong length
            batch_size=1,
        )


def test_compute_pareto_front():
    """Test Pareto front computation."""
    pytest.importorskip("torch")
    pytest.importorskip("botorch")
    
    from amprenta_rag.analysis.bayesian_optimization import compute_pareto_front
    
    # Simple 2D objectives where (1, 1) dominates all others
    objectives = [
        [0.5, 0.5],
        [0.7, 0.3],
        [0.3, 0.7],
        [1.0, 1.0],  # Pareto optimal
    ]
    
    result = compute_pareto_front(objectives)
    
    assert "indices" in result
    assert "values" in result
    # Index 3 should be on Pareto front (highest in both objectives)
    assert 3 in result["indices"]
    assert len(result["values"]) == len(result["indices"])


