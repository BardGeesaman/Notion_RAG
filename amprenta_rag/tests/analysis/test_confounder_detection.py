from __future__ import annotations

import pandas as pd
import pytest

from amprenta_rag.analysis import confounder_detection as cd


class _FakeStats:
    def chi2_contingency(self, contingency):
        return (1.0, 0.01, None, None)

    def ttest_ind(self, a, b):
        return (0.5, 0.02)

    def f_oneway(self, *groups):
        return (1.1, 0.03)


def test_detect_confounders_categorical(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", True)
    monkeypatch.setattr(cd, "stats", _FakeStats())

    df = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "cat": ["x", "y", "x", "y"],
            "num": [1, 2, 3, 4],
        }
    )

    results = cd.detect_confounders(df, "group")

    assert any(r["column"] == "cat" for r in results)
    assert any(r["column"] == "num" for r in results)


def test_detect_confounders_group_missing(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", True)
    monkeypatch.setattr(cd, "stats", _FakeStats())

    df = pd.DataFrame({"cat": ["a", "b"]})
    assert cd.detect_confounders(df, "group") == []


def test_detect_confounders_needs_scipy(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", False)
    with pytest.raises(ImportError):
        cd.detect_confounders(pd.DataFrame({"group": ["a", "b"]}), "group")


def test_detect_confounders_insufficient_groups(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", True)
    monkeypatch.setattr(cd, "stats", _FakeStats())
    df = pd.DataFrame({"group": ["A", "A"], "cat": ["x", "y"]})
    assert cd.detect_confounders(df, "group") == []


def test_get_confounder_report(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", True)
    monkeypatch.setattr(cd, "stats", _FakeStats())
    monkeypatch.setattr(cd, "calculate_imbalance_score", lambda series, groups: 0.4)

    df = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "cat": ["x", "y", "x", "y"],
            "num": [1, 2, 3, 4],
        }
    )

    report = cd.get_confounder_report(df, "group")

    assert report["confounders"]
    assert any("imbalance" in w.lower() for w in report["warnings"])
    assert "confounder" in report["summary"].lower()


def test_get_confounder_report_scipy_missing(monkeypatch) -> None:
    monkeypatch.setattr(cd, "SCIPY_AVAILABLE", False)
    report = cd.get_confounder_report(pd.DataFrame({"group": []}), "group")
    assert report["warnings"]
    assert "unavailable" in report["summary"].lower()


def test_calculate_imbalance_score_handles_small(monkeypatch) -> None:
    groups = pd.Series(["A"], index=[0])
    series = pd.Series(["x"], index=[0])
    assert cd.calculate_imbalance_score(series, groups) == 0.0


def test_calculate_imbalance_score_exception(monkeypatch) -> None:
    def _boom(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(cd.pd, "crosstab", _boom)
    groups = pd.Series(["A", "B"])
    series = pd.Series(["x", "y"])
    assert cd.calculate_imbalance_score(series, groups) == 0.0

