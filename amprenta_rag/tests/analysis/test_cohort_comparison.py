from __future__ import annotations

import numpy as np
import pandas as pd

from amprenta_rag.analysis import cohort_comparison as cc


def test_compute_effect_size_cases():
    assert np.isnan(cc.compute_effect_size([], []))
    assert cc.compute_effect_size([1, 1], [1, 1]) == 0.0
    val = cc.compute_effect_size([1, 2, 3], [4, 5, 6])
    assert val < 0  # second group larger


def test_compare_cohorts_uses_anova_and_adjust(monkeypatch):
    # Monkeypatch statistical functions for determinism
    monkeypatch.setattr(cc, "anova_oneway", lambda *groups: {"pvalue": 0.04, "statistic": 1.0})
    monkeypatch.setattr(cc, "adjust_pvalues", lambda pvals: (None, [p * 0.5 for p in pvals]))

    df = pd.DataFrame([[1, 2, 3], [3, 2, 1]], index=["f1", "f2"])
    cohort_labels = ["A", "B", "B"]
    out = cc.compare_cohorts(df, cohort_labels)
    assert set(out["feature"]) == {"f1", "f2"}
    assert (out["adj_pvalue"] <= 0.02).all()


def test_run_pairwise_comparisons(monkeypatch):
    monkeypatch.setattr(cc, "ttest_independent", lambda g1, g2: {"pvalue": 0.03})
    monkeypatch.setattr(cc, "adjust_pvalues", lambda pvals: (None, pvals))

    df = pd.DataFrame([[1, 2, 3, 4]], index=["feat"])
    labels = ["A", "A", "B", "B"]
    out = cc.run_pairwise_comparisons(df, labels)
    assert out.iloc[0]["feature"] == "feat"
    assert out.iloc[0]["interpretation"] in {"small", "moderate+", "large", "ns"}
    assert "cohens_d" in out
from __future__ import annotations

import numpy as np
import pandas as pd

from amprenta_rag.analysis import cohort_comparison as cc


def test_compute_effect_size_small_samples():
    d = cc.compute_effect_size([1, 2], [1, 2])
    assert d == 0 or np.isnan(d)


def test_compare_cohorts_adds_adj_pvalue():
    df = pd.DataFrame([[1, 2, 3, 4]], index=["feat"])
    labels = ["A", "A", "B", "B"]
    out = cc.compare_cohorts(df, labels)
    assert "adj_pvalue" in out.columns
    assert len(out) == 1


def test_run_pairwise_comparisons_returns_pairs():
    df = pd.DataFrame([[1, 2, 3, 4]], index=["feat"])
    labels = ["A", "A", "B", "B"]
    out = cc.run_pairwise_comparisons(df, labels)
    assert {"feature", "cohort1", "cohort2", "pvalue", "cohens_d", "interpretation", "adj_pvalue"}.issubset(out.columns)
    assert len(out) == 1

