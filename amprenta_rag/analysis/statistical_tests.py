from __future__ import annotations

from typing import Dict, Iterable, Sequence, Tuple

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _to_array(values: Iterable[float]) -> np.ndarray:
    arr = np.array(list(values), dtype=float)
    return arr[~np.isnan(arr)]


def _validate_pairs(a: np.ndarray, b: np.ndarray, min_len: int = 1) -> bool:
    return len(a) >= min_len and len(b) >= min_len


def ttest_independent(group1: Iterable[float], group2: Iterable[float]) -> Dict[str, float | str]:
    a = _to_array(group1)
    b = _to_array(group2)
    if not _validate_pairs(a, b, min_len=2):
        return {"statistic": np.nan, "pvalue": np.nan, "interpretation": "insufficient data"}
    stat, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
    return _format_result(stat, p)


def anova_oneway(*groups: Iterable[float]) -> Dict[str, float | str]:
    arrays = [g for g in (_to_array(g) for g in groups) if len(g) > 1]
    if len(arrays) < 2:
        return {"statistic": np.nan, "pvalue": np.nan, "interpretation": "insufficient data"}
    stat, p = stats.f_oneway(*arrays)
    return _format_result(stat, p)


def mann_whitney(group1: Iterable[float], group2: Iterable[float]) -> Dict[str, float | str]:
    a = _to_array(group1)
    b = _to_array(group2)
    if not _validate_pairs(a, b, min_len=1):
        return {"statistic": np.nan, "pvalue": np.nan, "interpretation": "insufficient data"}
    stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return _format_result(stat, p)


def pearson_corr(x: Iterable[float], y: Iterable[float]) -> Dict[str, float | str]:
    a = _to_array(x)
    b = _to_array(y)
    if len(a) < 2 or len(b) < 2 or len(a) != len(b):
        return {"statistic": np.nan, "pvalue": np.nan, "interpretation": "insufficient data"}
    r, p = stats.pearsonr(a, b)
    return _format_result(r, p)


def adjust_pvalues(pvalues: Sequence[float], method: str = "fdr_bh") -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply multiple testing correction.
    method: 'fdr_bh', 'bonferroni', etc. (statsmodels multipletests)
    Returns (adjusted_pvals, reject_flags)
    """
    arr = np.array(pvalues, dtype=float)
    valid = ~np.isnan(arr)
    adj = np.full_like(arr, np.nan)
    reject = np.zeros_like(arr, dtype=bool)
    if valid.sum() == 0:
        return adj, reject
    # multipletests returns (reject, pvals_corrected, alphacSidak, alphacBonf)
    reject_valid, adj_valid, _, _ = multipletests(arr[valid], method=method)
    adj[valid] = adj_valid
    reject[valid] = reject_valid
    return adj, reject


def _format_result(stat: float, p: float) -> Dict[str, float | str]:
    interp = "not significant"
    if not np.isnan(p):
        if p < 0.001:
            interp = "*** significant"
        elif p < 0.01:
            interp = "** significant"
        elif p < 0.05:
            interp = "* significant"
    return {"statistic": float(stat), "pvalue": float(p), "interpretation": interp}

