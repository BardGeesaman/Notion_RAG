"""Tests for compound ranking scorer module."""

from __future__ import annotations

import pytest

from amprenta_rag.ml.ranking.scorer import (
    ObjectiveScore,
    compute_pareto_ranks,
    compute_weighted_score,
    normalize_alerts,
    normalize_herg,
    normalize_logp,
    normalize_logs,
    normalize_potency,
)


def test_normalize_potency_nm() -> None:
    """Test potency normalization with nM units."""
    # 100 nM = pIC50 of 7 = normalized score of (7-3)/7 = 0.571
    assert abs(normalize_potency(100.0, "nM") - 0.571) < 0.01
    
    # 1 nM = pIC50 of 9 = normalized score of (9-3)/7 = 0.857
    assert abs(normalize_potency(1.0, "nM") - 0.857) < 0.01
    
    # Very high potency should cap at 1.0
    assert normalize_potency(0.01, "nM") >= 0.95


def test_normalize_potency_um() -> None:
    """Test potency normalization with µM units."""
    # 1 µM = 1000 nM = pIC50 of 6 = normalized score of (6-3)/7 = 0.429
    assert abs(normalize_potency(1.0, "µM") - 0.429) < 0.01
    assert abs(normalize_potency(1.0, "uM") - 0.429) < 0.01  # Alternative spelling
    
    # 10 µM = 10000 nM = pIC50 of 5 = normalized score of (5-3)/7 = 0.286
    assert abs(normalize_potency(10.0, "µM") - 0.286) < 0.01


def test_normalize_herg() -> None:
    """Test hERG normalization (inverted - lower blocking = higher safety)."""
    # Low hERG blocking = high safety score
    assert normalize_herg(0.1) == pytest.approx(0.9, rel=1e-9)
    assert normalize_herg(0.0) == pytest.approx(1.0, rel=1e-9)
    
    # High hERG blocking = low safety score
    assert normalize_herg(0.9) == pytest.approx(0.1, rel=1e-9)
    assert normalize_herg(1.0) == pytest.approx(0.0, rel=1e-9)
    
    # Values outside 0-1 should be clamped
    assert normalize_herg(-0.1) == pytest.approx(1.0, rel=1e-9)
    assert normalize_herg(1.5) == pytest.approx(0.0, rel=1e-9)


def test_normalize_logs() -> None:
    """Test logS normalization (higher = better solubility)."""
    # logS of -3 should give ((-3) + 6) / 6 = 0.5
    assert abs(normalize_logs(-3.0) - 0.5) < 0.01
    
    # logS of 0 (highly soluble) should give 1.0
    assert normalize_logs(0.0) == 1.0
    
    # logS of -6 (poorly soluble) should give 0.0
    assert normalize_logs(-6.0) == 0.0
    
    # Values outside range should be clamped
    assert normalize_logs(-10.0) == 0.0
    assert normalize_logs(2.0) == 1.0


def test_normalize_logp() -> None:
    """Test logP normalization (ideal range 1-3, target=2)."""
    # logP of 2 (ideal) should give 1.0
    assert normalize_logp(2.0) == 1.0
    
    # logP of 1 or 3 (edges of ideal range) should give 0.75
    assert abs(normalize_logp(1.0) - 0.75) < 0.01
    assert abs(normalize_logp(3.0) - 0.75) < 0.01
    
    # logP of 0 or 4 (deviation of 2) should give 0.5
    assert abs(normalize_logp(0.0) - 0.5) < 0.01
    assert abs(normalize_logp(4.0) - 0.5) < 0.01
    
    # Extreme values should not go negative
    assert normalize_logp(-5.0) >= 0.0
    assert normalize_logp(10.0) >= 0.0


def test_normalize_alerts_pains_weighted() -> None:
    """Test structural alerts normalization with PAINS weighting."""
    # No alerts = perfect safety score
    assert normalize_alerts([]) == 1.0
    
    # Single PAINS alert (weight 1.0) should give exp(-1/5) ≈ 0.819
    pains_alert = [{"type": "PAINS", "count": 1}]
    assert abs(normalize_alerts(pains_alert) - 0.819) < 0.01
    
    # Single Brenk alert (weight 0.5) should give exp(-0.5/5) ≈ 0.905
    brenk_alert = [{"type": "BRENK", "count": 1}]
    assert abs(normalize_alerts(brenk_alert) - 0.905) < 0.01
    
    # Multiple alerts should accumulate
    mixed_alerts = [
        {"type": "PAINS", "count": 1},
        {"type": "BRENK", "count": 2}
    ]
    # Total weight = 1.0 + 2*0.5 = 2.0, score = exp(-2/5) ≈ 0.670
    assert abs(normalize_alerts(mixed_alerts) - 0.670) < 0.01


def test_compute_weighted_score() -> None:
    """Test weighted score computation."""
    objectives = [
        ObjectiveScore("potency", 100.0, 0.8, 0.4),
        ObjectiveScore("herg", 0.2, 0.8, 0.3),
        ObjectiveScore("alerts", 0, 1.0, 0.3)
    ]
    
    # Expected: (0.8*0.4 + 0.8*0.3 + 1.0*0.3) / (0.4+0.3+0.3) = 0.86
    expected = (0.8 * 0.4 + 0.8 * 0.3 + 1.0 * 0.3) / 1.0
    assert abs(compute_weighted_score(objectives) - expected) < 0.01


def test_compute_pareto_ranks_simple() -> None:
    """Test Pareto ranking with simple dominance."""
    compounds = [
        {"id": "A", "objectives": {"potency": 0.9, "safety": 0.8}},  # Dominates B
        {"id": "B", "objectives": {"potency": 0.7, "safety": 0.6}},  # Dominated by A
        {"id": "C", "objectives": {"potency": 0.6, "safety": 0.9}},  # Non-dominated
    ]
    
    ranks = compute_pareto_ranks(compounds)
    
    # A and C should be rank 1 (Pareto front), B should be rank 2
    assert ranks[0] == 1  # A
    assert ranks[1] == 2  # B (dominated by A)
    assert ranks[2] == 1  # C (non-dominated)


def test_compute_pareto_ranks_dominated() -> None:
    """Test Pareto ranking with clear dominance hierarchy."""
    compounds = [
        {"id": "Best", "objectives": {"obj1": 1.0, "obj2": 1.0}},    # Rank 1
        {"id": "Good", "objectives": {"obj1": 0.8, "obj2": 0.8}},    # Rank 2
        {"id": "Poor", "objectives": {"obj1": 0.5, "obj2": 0.5}},    # Rank 3
    ]
    
    ranks = compute_pareto_ranks(compounds)
    
    assert ranks[0] == 1  # Best
    assert ranks[1] == 2  # Good
    assert ranks[2] == 3  # Poor


def test_compute_pareto_ranks_empty() -> None:
    """Test Pareto ranking with empty input."""
    assert compute_pareto_ranks([]) == []


def test_compute_pareto_ranks_single() -> None:
    """Test Pareto ranking with single compound."""
    compounds = [{"id": "Only", "objectives": {"obj1": 0.5}}]
    ranks = compute_pareto_ranks(compounds)
    assert ranks == [1]


def test_normalize_potency_edge_cases() -> None:
    """Test potency normalization edge cases."""
    # Zero or negative values should return 0
    assert normalize_potency(0.0, "nM") == 0.0
    assert normalize_potency(-1.0, "nM") == 0.0
    
    # Unknown units should default to nM with warning
    assert normalize_potency(100.0, "unknown") > 0.0


def test_normalize_alerts_custom_weights() -> None:
    """Test structural alerts with custom weights."""
    alerts = [{"type": "CUSTOM", "count": 1}]
    custom_weights = {"CUSTOM": 2.0}
    
    # Should use custom weight of 2.0
    score = normalize_alerts(alerts, custom_weights)
    expected = pytest.approx(0.670, abs=0.01)  # exp(-2/5)
    assert score == expected


def test_compute_weighted_score_empty() -> None:
    """Test weighted score with no objectives."""
    assert compute_weighted_score([]) == 0.0


def test_compute_weighted_score_zero_weights() -> None:
    """Test weighted score with zero total weight."""
    objectives = [
        ObjectiveScore("test", 1.0, 0.5, 0.0)
    ]
    assert compute_weighted_score(objectives) == 0.0
