"""Confounder detection utilities for experimental design analysis."""
from __future__ import annotations

from typing import List, Dict, Any

import pandas as pd

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    stats = None  # type: ignore


def detect_confounders(
    metadata_df: pd.DataFrame,
    group_column: str,
) -> List[Dict[str, Any]]:
    """
    Detect potential confounders in metadata by testing associations with group variable.

    Args:
        metadata_df: DataFrame with metadata columns
        group_column: Name of column containing group labels

    Returns:
        List of dicts with column, test, p_value, is_confounder
    """
    if not SCIPY_AVAILABLE:
        logger.error("[CONFOUNDER] scipy.stats not available")
        raise ImportError("scipy.stats is required for confounder detection")

    if group_column not in metadata_df.columns:
        logger.error("[CONFOUNDER] Group column '%s' not found", group_column)
        return []

    results = []
    groups = metadata_df[group_column].dropna()

    if len(groups.unique()) < 2:
        logger.warning("[CONFOUNDER] Need at least 2 groups for confounder detection")
        return []

    for col in metadata_df.columns:
        if col == group_column:
            continue

        series = metadata_df[col].dropna()
        if len(series) == 0:
            continue

        # Align series with groups
        aligned_idx = series.index.intersection(groups.index)
        if len(aligned_idx) < 2:
            continue

        aligned_series = series.loc[aligned_idx]
        aligned_groups = groups.loc[aligned_idx]

        try:
            # Determine if categorical or numeric
            is_categorical = (
                aligned_series.dtype == 'object' or
                aligned_series.dtype.name == 'category' or
                aligned_series.nunique() < 10  # Heuristic: few unique values = categorical
            )

            if is_categorical:
                # Chi-square test
                try:
                    contingency = pd.crosstab(aligned_series, aligned_groups)
                    if contingency.size == 0:
                        continue

                    chi2, p_value = stats.chi2_contingency(contingency)[:2]
                    test_name = "chi-square"

                    logger.debug("[CONFOUNDER] %s vs %s: chi2=%.3f, p=%.4f", col, group_column, chi2, p_value)
                except Exception as e:
                    logger.warning("[CONFOUNDER] Chi-square failed for %s: %r", col, e)
                    continue

            else:
                # Numeric: use t-test or ANOVA
                group_values = {}
                for group in aligned_groups.unique():
                    group_values[group] = aligned_series[aligned_groups == group].dropna()

                if len(group_values) == 2:
                    # Two groups: t-test
                    groups_list = list(group_values.values())
                    if len(groups_list[0]) < 2 or len(groups_list[1]) < 2:
                        continue

                    t_stat, p_value = stats.ttest_ind(groups_list[0], groups_list[1])
                    test_name = "t-test"

                    logger.debug("[CONFOUNDER] %s vs %s: t=%.3f, p=%.4f", col, group_column, t_stat, p_value)

                else:
                    # Multiple groups: ANOVA
                    groups_list = [v.values for v in group_values.values() if len(v) >= 2]
                    if len(groups_list) < 2:
                        continue

                    f_stat, p_value = stats.f_oneway(*groups_list)
                    test_name = "ANOVA"

                    logger.debug("[CONFOUNDER] %s vs %s: F=%.3f, p=%.4f", col, group_column, f_stat, p_value)

            is_confounder = p_value < 0.05 if p_value is not None else False

            results.append({
                "column": col,
                "test": test_name,
                "p_value": float(p_value) if p_value is not None else None,
                "is_confounder": is_confounder,
            })

        except Exception as e:
            logger.warning("[CONFOUNDER] Error testing %s: %r", col, e)
            continue

    logger.info("[CONFOUNDER] Detected %d potential confounders out of %d variables",
               sum(r["is_confounder"] for r in results), len(results))
    return results


def calculate_imbalance_score(series: pd.Series, groups: pd.Series) -> float:
    """
    Calculate imbalance score (variance of group proportions).

    Higher score = more imbalanced distribution across groups.

    Args:
        series: Categorical series
        groups: Group labels

    Returns:
        Imbalance score (0.0 = balanced, higher = more imbalanced)
    """
    try:
        # Align series with groups
        aligned_idx = series.index.intersection(groups.index)
        if len(aligned_idx) < 2:
            return 0.0

        aligned_series = series.loc[aligned_idx]
        aligned_groups = groups.loc[aligned_idx]

        # Calculate proportions for each category within each group
        contingency = pd.crosstab(aligned_series, aligned_groups)
        proportions = contingency.div(contingency.sum(axis=0), axis=1)

        # Calculate variance of proportions across groups for each category
        variances = proportions.var(axis=1)

        # Average variance across categories
        imbalance_score = float(variances.mean()) if len(variances) > 0 else 0.0

        logger.debug("[CONFOUNDER] Imbalance score: %.4f", imbalance_score)
        return imbalance_score

    except Exception as e:
        logger.warning("[CONFOUNDER] Error calculating imbalance score: %r", e)
        return 0.0


def get_confounder_report(
    metadata_df: pd.DataFrame,
    group_column: str,
) -> Dict[str, Any]:
    """
    Generate a comprehensive confounder detection report.

    Args:
        metadata_df: DataFrame with metadata columns
        group_column: Name of column containing group labels

    Returns:
        Dict with confounders, warnings, and summary
    """
    if not SCIPY_AVAILABLE:
        return {
            "confounders": [],
            "warnings": ["scipy.stats not available"],
            "summary": "Confounder detection unavailable",
        }

    try:
        logger.info("[CONFOUNDER] Generating confounder report for group: %s", group_column)

        # Detect confounders
        confounder_results = detect_confounders(metadata_df, group_column)

        # Extract confounders
        confounders = [
            {
                "column": r["column"],
                "test": r["test"],
                "p_value": r["p_value"],
            }
            for r in confounder_results
            if r["is_confounder"]
        ]

        # Generate warnings
        warnings = []

        if len(confounders) > 0:
            warnings.append(f"Found {len(confounders)} potential confounder(s) (p < 0.05)")

        # Check for imbalance in categorical variables
        groups = metadata_df[group_column].dropna()
        for col in metadata_df.columns:
            if col == group_column:
                continue

            series = metadata_df[col].dropna()
            if series.dtype == 'object' or series.nunique() < 10:
                imbalance = calculate_imbalance_score(series, groups)
                if imbalance > 0.3:  # Threshold for high imbalance
                    warnings.append(f"High imbalance detected in '{col}' (score: {imbalance:.2f})")

        # Generate summary
        total_vars = len([c for c in metadata_df.columns if c != group_column])
        confounder_count = len(confounders)

        if confounder_count == 0:
            summary = f"✓ No confounders detected among {total_vars} variables"
        else:
            summary = f"⚠ {confounder_count} potential confounder(s) detected among {total_vars} variables"

        report = {
            "confounders": confounders,
            "warnings": warnings,
            "summary": summary,
        }

        logger.info("[CONFOUNDER] Report: %s", summary)
        return report

    except Exception as e:
        logger.error("[CONFOUNDER] Error generating report: %r", e)
        return {
            "confounders": [],
            "warnings": [f"Error: {str(e)}"],
            "summary": "Confounder detection failed",
        }
