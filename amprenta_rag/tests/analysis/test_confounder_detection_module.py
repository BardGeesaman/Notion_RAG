from __future__ import annotations

from typing import Any, Dict

import pandas as pd
import pytest

from amprenta_rag.analysis import confounder_detection as cd


def test_detect_confounders_missing_group_column_logs_and_returns_empty(caplog):
    df = pd.DataFrame({"a": [1, 2]})
    res = cd.detect_confounders(df, "group")
    assert res == []


def test_detect_confounders_categorical_and_numeric(monkeypatch):
    # Ensure scipy available path; mock stats
    class FakeStats:
        @staticmethod
        def chi2_contingency(cont):
            return (1.0, 0.04, None, None)

        @staticmethod
        def ttest_ind(a, b):
            return (0.0, 0.03)

        @staticmethod
        def f_oneway(*args):
            return (1.0, 0.02)

    monkeypatch.setitem(cd.__dict__, "stats", FakeStats())
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    df = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "cat": ["x", "x", "y", "y"],
            "num": [1, 2, 3, 4],
        }
    )
    res = cd.detect_confounders(df, "group")
    fields = {r["column"] for r in res}
    assert {"cat", "num"}.issubset(fields)
    assert any(r["is_confounder"] for r in res)


def test_detect_confounders_insufficient_groups(caplog, monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)
    df = pd.DataFrame({"group": ["A", "A"], "x": [1, 2]})
    res = cd.detect_confounders(df, "group")
    assert res == []


def test_calculate_imbalance_score_handles_small_and_errors(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "logger", type("L", (), {"warning": lambda *a, **k: None})())
    series = pd.Series(["x", "y"])
    groups = pd.Series(["A", "A"])
    score = cd.calculate_imbalance_score(series, groups)
    assert score == 0.0


def test_get_confounder_report_with_scipy(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    def fake_detect(df, grp):
        return [
            {"column": "c1", "test": "chi", "p_value": 0.01, "is_confounder": True},
            {"column": "c2", "test": "ttest", "p_value": 0.2, "is_confounder": False},
        ]

    monkeypatch.setitem(cd.__dict__, "detect_confounders", fake_detect)

    df = pd.DataFrame({"group": ["A", "B"], "c1": [1, 2]})
    report: Dict[str, Any] = cd.get_confounder_report(df, "group")
    assert report["confounders"][0]["column"] == "c1"
    assert "warnings" in report


def test_get_confounder_report_without_scipy(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", False)
    report = cd.get_confounder_report(pd.DataFrame(), "group")
    assert report["warnings"]


def test_detect_confounders_raises_without_scipy(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", False)
    with pytest.raises(ImportError):
        cd.detect_confounders(pd.DataFrame(), "group")


def test_get_confounder_report_imbalance_warning(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    df = pd.DataFrame(
        {
            "group": ["A"] * 5 + ["B"] * 5,
            "cat": ["x"] * 5 + ["y"] * 5,
        }
    )

    # reuse real detect_confounders but speed by bypassing stats functions
    class FakeStats:
        @staticmethod
        def chi2_contingency(cont):
            return (1.0, 0.1, None, None)

        @staticmethod
        def ttest_ind(a, b):
            return (0.0, 0.5)

        @staticmethod
        def f_oneway(*args):
            return (0.0, 0.5)

    monkeypatch.setitem(cd.__dict__, "stats", FakeStats())
    report = cd.get_confounder_report(df, "group")
    assert any("imbalance" in w for w in report["warnings"])


def test_detect_confounders_handles_chi_square_failure(monkeypatch):
    class FakeStats:
        @staticmethod
        def chi2_contingency(cont):
            raise RuntimeError("bad chi")

        @staticmethod
        def ttest_ind(a, b):
            return (0.0, 0.5)

        @staticmethod
        def f_oneway(*args):
            return (0.0, 0.5)

    monkeypatch.setitem(cd.__dict__, "stats", FakeStats())
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)
    df = pd.DataFrame({"group": ["A", "B"], "cat": ["x", "y"]})
    res = cd.detect_confounders(df, "group")
    assert res == []


def test_get_confounder_report_no_confounders_summary(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)
    monkeypatch.setitem(cd.__dict__, "detect_confounders", lambda df, grp: [])
    df = pd.DataFrame({"group": ["A", "B"], "x": [1, 2]})
    report = cd.get_confounder_report(df, "group")
    assert "No confounders" in report["summary"]


def test_get_confounder_report_error_path(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    def boom(*a, **k):
        raise RuntimeError("fail")

    monkeypatch.setitem(cd.__dict__, "detect_confounders", boom)
    df = pd.DataFrame({"group": ["A", "B"], "x": [1, 2]})
    report = cd.get_confounder_report(df, "group")
    assert "Error:" in report["warnings"][0]


def test_detect_confounders_anova_path(monkeypatch):
    class FakeStats:
        @staticmethod
        def chi2_contingency(cont):
            return (1.0, 0.6, None, None)

        @staticmethod
        def ttest_ind(a, b):
            return (0.0, 0.6)

        @staticmethod
        def f_oneway(*args):
            return (2.0, 0.01)

    monkeypatch.setitem(cd.__dict__, "stats", FakeStats())
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    df = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B", "C", "C"],
            "num": [1, 2, 3, 4, 5, 6],
        }
    )
    res = cd.detect_confounders(df, "group")
    assert any(r["test"] == "ANOVA" and r["is_confounder"] for r in res)


def test_get_confounder_report_with_confounders(monkeypatch):
    monkeypatch.setitem(cd.__dict__, "SCIPY_AVAILABLE", True)

    def fake_detect(df, grp):
        return [
            {"column": "c1", "test": "chi", "p_value": 0.01, "is_confounder": True},
            {"column": "c2", "test": "ttest", "p_value": 0.2, "is_confounder": False},
        ]

    monkeypatch.setitem(cd.__dict__, "detect_confounders", fake_detect)
    df = pd.DataFrame({"group": ["A", "B"], "c1": [1, 2], "c2": [0, 1]})
    report = cd.get_confounder_report(df, "group")
    assert "potential confounder" in report["summary"]
    assert any("Found 1" in w for w in report["warnings"])

