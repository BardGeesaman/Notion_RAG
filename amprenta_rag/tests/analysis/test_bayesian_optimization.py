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


