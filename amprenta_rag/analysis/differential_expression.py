from __future__ import annotations

from typing import Iterable, List

import numpy as np
import pandas as pd

from amprenta_rag.analysis.statistical_tests import adjust_pvalues, ttest_independent


def classify_significance(row: pd.Series, fc_threshold: float, pval_threshold: float) -> str:
    if pd.isna(row["log2_fc"]) or pd.isna(row["adj_pvalue"]):
        return "ns"
    if row["adj_pvalue"] > pval_threshold:
        return "ns"
    if row["log2_fc"] >= fc_threshold:
        return "up"
    if row["log2_fc"] <= -fc_threshold:
        return "down"
    return "ns"


def run_differential_expression(
    df: pd.DataFrame,
    group1_idx: Iterable[int],
    group2_idx: Iterable[int],
    fc_threshold: float = 1.5,
    pval_threshold: float = 0.05,
) -> pd.DataFrame:
    """Compute differential expression between two groups.

    Args:
        df: DataFrame of shape (features x samples) with numeric values.
        group1_idx: indices of samples in group 1
        group2_idx: indices of samples in group 2
        fc_threshold: absolute log2 fold-change threshold for significance
        pval_threshold: adjusted p-value threshold

    Returns:
        DataFrame with columns: feature, log2_fc, pvalue, adj_pvalue, significant
    """
    features: List[str] = df.index.tolist()
    data = df.to_numpy(dtype=float)

    g1 = list(group1_idx)
    g2 = list(group2_idx)

    log2_fc_list: List[float] = []
    pvals: List[float] = []

    for row in data:
        g1_vals = row[g1]
        g2_vals = row[g2]
        # fold change: mean difference
        mean1 = np.nanmean(g1_vals)
        mean2 = np.nanmean(g2_vals)
        if mean2 == 0 or np.isnan(mean1) or np.isnan(mean2):
            log2_fc = np.nan
        else:
            fc = (mean1 + 1e-9) / (mean2 + 1e-9)
            log2_fc = np.log2(fc)
        log2_fc_list.append(log2_fc)

        res = ttest_independent(g1_vals, g2_vals)
        pval = res.get("pvalue", np.nan)
        pvals.append(float(pval) if pval is not None else np.nan)

    adj_reject, adj_pvals = adjust_pvalues(pvals)

    out = pd.DataFrame(
        {
            "feature": features,
            "log2_fc": log2_fc_list,
            "pvalue": pvals,
            "adj_pvalue": adj_pvals,
        }
    )
    out["significant"] = out.apply(lambda r: classify_significance(r, fc_threshold, pval_threshold), axis=1)
    return out


__all__ = ["run_differential_expression", "classify_significance"]
