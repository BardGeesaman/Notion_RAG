from __future__ import annotations

import pytest

from amprenta_rag.analysis import power_analysis as pa


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

