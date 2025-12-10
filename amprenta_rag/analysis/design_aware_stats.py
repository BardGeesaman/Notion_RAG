from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd

from amprenta_rag.analysis.statistical_tests import (
    adjust_pvalues,
    anova_oneway,
    pearson_corr,
    ttest_independent,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _extract_group_values(features_df: pd.DataFrame, sample_ids: List[str]) -> np.ndarray:
    existing = [s for s in sample_ids if s in features_df.columns]
    if not existing:
        return np.array([])
    return features_df[existing].to_numpy()


def _ttest_by_feature(features_df: pd.DataFrame, group_a: List[str], group_b: List[str]) -> Dict[str, float]:
    pvals: Dict[str, float] = {}
    values_a = _extract_group_values(features_df, group_a)
    values_b = _extract_group_values(features_df, group_b)
    if values_a.size == 0 or values_b.size == 0:
        return {f: np.nan for f in features_df.index}
    for idx, feature in enumerate(features_df.index):
        res = ttest_independent(values_a[idx, :], values_b[idx, :])
        pvals[feature] = res.get("pvalue", np.nan)
    return pvals


def _anova_by_feature(features_df: pd.DataFrame, groups: Dict[str, List[str]]) -> Dict[str, float]:
    pvals: Dict[str, float] = {}
    # Pre-extract arrays per group
    group_arrays: Dict[str, np.ndarray] = {}
    for name, ids in groups.items():
        arr = _extract_group_values(features_df, ids)
        if arr.size > 0:
            group_arrays[name] = arr
    for idx, feature in enumerate(features_df.index):
        arrays = [arr[idx, :] for arr in group_arrays.values() if arr.shape[0] > idx]
        res = anova_oneway(*arrays)
        pvals[feature] = res.get("pvalue", np.nan)
    return pvals


def _dose_response_by_feature(features_df: pd.DataFrame, groups: Dict[str, List[str]]) -> Dict[str, float]:
    """
    Correlate feature means against dose ordering (best-effort parsing of group names).
    """
    pvals: Dict[str, float] = {}
    # Determine dose ordering
    dose_levels: List[str] = list(groups.keys())
    parsed = []
    for g in dose_levels:
        match = None
        try:
            match = float(g)
        except Exception:
            m = None
            try:
                import re as _re
                m = _re.search(r"([\d.]+)", g)
            except Exception:
                m = None
            if m:
                try:
                    match = float(m.group(1))
                except Exception:
                    match = None
        parsed.append(match)
    # Fallback to ordinal if parsing fails
    dose_values = []
    next_ord = 1.0
    for p in parsed:
        if p is None:
            dose_values.append(next_ord)
            next_ord += 1.0
        else:
            dose_values.append(p)

    for idx, feature in enumerate(features_df.index):
        means: List[float] = []
        doses: List[float] = []
        for dose, sample_ids in zip(dose_values, dose_levels):
            arr = _extract_group_values(features_df, groups[sample_ids])
            if arr.size == 0 or arr.shape[0] <= idx:
                continue
            vals = arr[idx, :]
            if vals.size == 0:
                continue
            means.append(float(np.nanmean(vals)))
            doses.append(dose)
        if len(means) >= 2 and len(doses) == len(means):
            res = pearson_corr(doses, means)
            pvals[feature] = res.get("pvalue", np.nan)
        else:
            pvals[feature] = np.nan
    return pvals


def run_design_aware_analysis(
    features_df: pd.DataFrame,
    sample_groups: Dict[str, List[str]],
    design_type: str,
) -> Dict[str, any]:
    """
    Run statistical analysis tailored to experimental design.

    Returns:
        dict with keys: test_type, p_values (raw), p_values_adj, significant_features, design_type
    """
    design = (design_type or "observational") if isinstance(design_type, str) else "observational"
    design = design.lower()
    test_type = "descriptive"
    p_values: Dict[str, float] = {}

    if design == "case_control":
        control = sample_groups.get("control") or sample_groups.get("controls") or []
        case = sample_groups.get("case") or sample_groups.get("cases") or []
        if not control and not case and len(sample_groups) >= 2:
            # Fallback: first two groups
            keys = list(sample_groups.keys())
            control, case = sample_groups.get(keys[0], []), sample_groups.get(keys[1], [])
        p_values = _ttest_by_feature(features_df, control, case)
        test_type = "t-test"

    elif design == "time_course":
        # Expect sample_groups like {"timepoints": {...}} or keyed by time labels
        if "timepoints" in sample_groups and isinstance(sample_groups["timepoints"], dict):
            groups = sample_groups["timepoints"]
        else:
            groups = sample_groups
        p_values = _anova_by_feature(features_df, groups)
        test_type = "anova"

    elif design == "intervention":
        treated = sample_groups.get("treated") or []
        untreated = sample_groups.get("untreated") or sample_groups.get("control") or []
        p_values = _ttest_by_feature(features_df, treated, untreated)
        test_type = "t-test"

    elif design == "dose_response":
        p_values = _dose_response_by_feature(features_df, sample_groups)
        test_type = "correlation"

    elif design == "multi_factorial":
        # Placeholder: two-way ANOVA not available without statsmodels design matrix
        logger.info("[STATS] Two-way ANOVA not implemented; returning NaN p-values")
        p_values = {f: np.nan for f in features_df.index}
        test_type = "two-way anova (not implemented)"

    else:
        # Observational / fallback
        logger.info("[STATS] Observational design: no hypothesis test performed")
        p_values = {f: np.nan for f in features_df.index}
        test_type = "none"

    # Adjust p-values
    raw_pvals = list(p_values.values())
    adj, reject = adjust_pvalues(raw_pvals, method="fdr_bh")
    p_values_adj = dict(zip(p_values.keys(), adj))
    significant_features = [feat for feat, rej in zip(p_values.keys(), reject) if rej]

    logger.info(
        "[STATS] design=%s test=%s features=%d significant=%d",
        design,
        test_type,
        len(p_values),
        len(significant_features),
    )

    return {
        "design_type": design,
        "test_type": test_type,
        "p_values": p_values,
        "p_values_adj": p_values_adj,
        "significant_features": significant_features,
    }

