from __future__ import annotations

from itertools import combinations
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from amprenta_rag.analysis.statistical_tests import (
    adjust_pvalues,
    anova_oneway,
    ttest_independent,
)


def compute_effect_size(group1: Iterable[float], group2: Iterable[float]) -> float:
    a = np.array(list(group1), dtype=float)
    b = np.array(list(group2), dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled_std = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    if pooled_std == 0:
        return 0.0
    return (a.mean() - b.mean()) / pooled_std


def compare_cohorts(df: pd.DataFrame, cohort_labels: List[str]) -> pd.DataFrame:
    """Run ANOVA across cohorts for each feature."""
    features = df.index.tolist()
    results: List[Dict[str, float | str]] = []

    # Group columns by cohort label
    cohorts = {}
    for idx, label in enumerate(cohort_labels):
        cohorts.setdefault(label, []).append(idx)

    for feat, row in df.iterrows():
        groups = [row.values[idxs] for idxs in cohorts.values()]
        res = anova_oneway(*groups)
        results.append({
            "feature": feat,
            "pvalue": res.get("pvalue", np.nan),
            "statistic": res.get("statistic", np.nan),
        })

    out = pd.DataFrame(results)
    if not out.empty:
        _, adj = adjust_pvalues(out["pvalue"].tolist())
        out["adj_pvalue"] = adj
    return out


def run_pairwise_comparisons(df: pd.DataFrame, cohort_labels: List[str]) -> pd.DataFrame:
    """All pairwise t-tests between cohorts with effect sizes."""
    cohorts = {}
    for idx, label in enumerate(cohort_labels):
        cohorts.setdefault(label, []).append(idx)

    pairs = list(combinations(cohorts.keys(), 2))
    records: List[Dict[str, float | str]] = []

    for feat, row in df.iterrows():
        for c1, c2 in pairs:
            g1 = row.values[cohorts[c1]]
            g2 = row.values[cohorts[c2]]
            res = ttest_independent(g1, g2)
            d = compute_effect_size(g1, g2)
            interpretation = "ns"
            if not np.isnan(res.get("pvalue", np.nan)) and res["pvalue"] < 0.05:
                if d >= 0.5:
                    interpretation = "moderate+" if d < 0.8 else "large"
                elif d <= -0.5:
                    interpretation = "moderate+" if d > -0.8 else "large"
                else:
                    interpretation = "small"
            records.append({
                "feature": feat,
                "cohort1": c1,
                "cohort2": c2,
                "pvalue": res.get("pvalue", np.nan),
                "cohens_d": d,
                "interpretation": interpretation,
            })

    out = pd.DataFrame(records)
    if not out.empty:
        _, adj = adjust_pvalues(out["pvalue"].tolist())
        out["adj_pvalue"] = adj
    return out


__all__ = [
    "compute_effect_size",
    "compare_cohorts",
    "run_pairwise_comparisons",
]
