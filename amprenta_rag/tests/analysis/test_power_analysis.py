from __future__ import annotations

import pytest

from amprenta_rag.analysis import power_analysis as pa
from amprenta_rag.analysis.power_analysis import (
    estimate_effect_size_from_data,
    calculate_plate_layout,
    estimate_experiment_cost,
)


class DummyPower:
    def __init__(self, result=10, power=0.9):
        self._result = result
        self._power = power

    def solve_power(self, **kwargs):
        return self._result

    def power(self, **kwargs):
        return self._power


def test_get_effect_size_preset():
    assert pa.get_effect_size_preset("small") == 0.2
    assert pa.get_effect_size_preset("LARGE") == 0.8
    assert pa.get_effect_size_preset("unknown") == 0.5


def test_calculate_sample_size_with_stubs(monkeypatch):
    monkeypatch.setattr(pa, "STATSMODELS_AVAILABLE", True)
    monkeypatch.setattr(pa, "TTestIndPower", lambda: DummyPower(result=12))
    monkeypatch.setattr(pa, "FTestAnovaPower", lambda: DummyPower(result=30))
    monkeypatch.setattr(pa, "TTestPower", lambda: DummyPower(result=40))

    assert pa.calculate_sample_size(effect_size=0.5, test_type="t-test") == 12
    assert pa.calculate_sample_size(effect_size=0.5, test_type="anova") == 30
    assert pa.calculate_sample_size(effect_size=0.3, test_type="correlation") == 40
    assert pa.calculate_sample_size(effect_size=0.3, test_type="chi-square") == 40
    with pytest.raises(ValueError):
        pa.calculate_sample_size(effect_size=0.5, test_type="invalid")  # type: ignore[arg-type]


def test_calculate_power_with_stubs(monkeypatch):
    monkeypatch.setattr(pa, "STATSMODELS_AVAILABLE", True)
    monkeypatch.setattr(pa, "TTestIndPower", lambda: DummyPower(power=0.8))
    monkeypatch.setattr(pa, "FTestAnovaPower", lambda: DummyPower(power=0.7))
    monkeypatch.setattr(pa, "TTestPower", lambda: DummyPower(power=0.6))

    assert pa.calculate_power(n=10, effect_size=0.5, test_type="t-test") == 0.8
    assert pa.calculate_power(n=12, effect_size=0.4, test_type="anova") == 0.7
    assert pa.calculate_power(n=15, effect_size=0.3, test_type="correlation") == 0.6
    assert pa.calculate_power(n=15, effect_size=0.3, test_type="chi-square") == 0.6
    with pytest.raises(ValueError):
        pa.calculate_power(n=10, effect_size=0.5, test_type="invalid")  # type: ignore[arg-type]


def test_estimate_effect_size_from_data():
    """Test effect size estimation from data."""
    group1 = [10, 12, 11, 13, 12]
    group2 = [8, 9, 7, 8, 9]
    
    effect_size = estimate_effect_size_from_data(group1, group2)
    
    assert effect_size > 0, "Effect size should be positive"
    assert abs(effect_size) < 5, "Effect size should be reasonable"


def test_estimate_effect_size_empty_data():
    """Test effect size with empty data."""
    result = estimate_effect_size_from_data([], [1, 2, 3])
    assert result == 0.0
    
    result2 = estimate_effect_size_from_data([1, 2, 3], [])
    assert result2 == 0.0


def test_calculate_plate_layout_96():
    """Test plate layout for 96-well format."""
    result = calculate_plate_layout(100, plate_format=96)
    
    assert result["plates_needed"] == 2
    assert result["wells_used"] == 100
    assert result["empty_wells"] == 92


def test_calculate_plate_layout_384():
    """Test plate layout for 384-well format."""
    result = calculate_plate_layout(400, plate_format=384)
    
    assert result["plates_needed"] == 2
    assert result["wells_used"] == 400
    assert result["empty_wells"] == 368


def test_estimate_experiment_cost():
    """Test experiment cost estimation."""
    result = estimate_experiment_cost(n=100, cost_per_sample=10.0, overhead_pct=0.15)
    
    assert result["sample_cost"] == 1000.0
    assert result["overhead"] == 150.0
    assert result["total"] == 1150.0

