from __future__ import annotations

import pandas as pd

from amprenta_rag.analysis import differential_expression as de


def test_classify_significance():
    row = pd.Series({"log2_fc": 2.0, "adj_pvalue": 0.01})
    assert de.classify_significance(row, fc_threshold=1.0, pval_threshold=0.05) == "up"
    row = pd.Series({"log2_fc": -2.0, "adj_pvalue": 0.01})
    assert de.classify_significance(row, fc_threshold=1.0, pval_threshold=0.05) == "down"
    row = pd.Series({"log2_fc": 0.5, "adj_pvalue": 0.2})
    assert de.classify_significance(row, fc_threshold=1.0, pval_threshold=0.05) == "ns"


def test_run_differential_expression_basic():
    df = pd.DataFrame(
        [[1, 2, 10, 12], [5, 5, 5, 5]],
        index=["f1", "f2"],
    )
    res = de.run_differential_expression(df, group1_idx=[0, 1], group2_idx=[2, 3], fc_threshold=0.5, pval_threshold=0.1)
    assert set(res.columns) == {"feature", "log2_fc", "pvalue", "adj_pvalue", "significant"}
    assert len(res) == 2
    assert set(res["significant"].unique()).issubset({"up", "down", "ns"})

