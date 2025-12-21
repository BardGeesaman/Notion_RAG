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

