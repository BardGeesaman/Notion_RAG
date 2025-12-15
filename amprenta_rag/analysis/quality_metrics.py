from __future__ import annotations

from typing import Dict, List

from amprenta_rag.database.models import Dataset, Feature
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _feature_stats(feature: Feature) -> Dict[str, float]:
    """Extract fold-change and p-value stats from feature.external_ids."""
    data = feature.external_ids or {}
    if not isinstance(data, dict):
        return {}

    candidates = [data]
    for key in ("stats", "de_stats", "diffexp", "differential_expression"):
        nested = data.get(key)
        if isinstance(nested, dict):
            candidates.append(nested)

    def first(keys):
        for obj in candidates:
            for k in keys:
                if k in obj and obj[k] is not None:
                    try:
                        return float(obj[k])
                    except Exception:
                        continue
        return None

    fc = first(["log2FC", "log2fc", "log2_fold_change", "fold_change", "fc"])
    p = first(["pvalue", "p_value", "adj_pvalue", "adj_p", "pval", "p_val", "padj"])
    return {"log2FC": fc, "pvalue": p}


def compute_quality_score(dataset: Dataset) -> Dict[str, object]:
    """
    Compute a quality score for a dataset using simple heuristics.

    Returns:
        dict with keys: score (0-100), status ("high"/"medium"/"low"),
        issues (list[str]), metrics (dict).
    """
    issues: List[str] = []
    score = 100.0

    features = list(dataset.features or [])
    feature_count = len(features)
    if feature_count == 0:
        issues.append("No features linked to dataset")
        score -= 50

    with_stats = 0
    outliers = 0
    for feat in features:
        stats = _feature_stats(feat)
        if stats.get("log2FC") is not None and stats.get("pvalue") is not None:
            with_stats += 1
            try:
                if abs(stats["log2FC"]) > 5:
                    outliers += 1
            except Exception as e:
                logger.warning(
                    "[QUALITY] Skipping outlier check for feature %s: %r",
                    getattr(feat, "id", "unknown"),
                    e,
                )

    stats_coverage = (with_stats / feature_count * 100) if feature_count else 0
    if stats_coverage < 30:
        issues.append("Low stats coverage (<30%)")
        score -= 25
    elif stats_coverage < 70:
        issues.append("Moderate stats coverage (<70%)")
        score -= 10

    if outliers > max(1, 0.05 * feature_count):
        issues.append("High outlier rate in log2FC")
        score -= 10

    # Missing metadata penalties
    if not dataset.description:
        issues.append("Missing description")
        score -= 5
    if not dataset.disease:
        issues.append("Missing disease annotation")
        score -= 5

    score = max(0.0, min(100.0, score))
    if score >= 80:
        status = "high"
    elif score >= 50:
        status = "medium"
    else:
        status = "low"

    metrics = {
        "feature_count": feature_count,
        "stats_coverage_pct": round(stats_coverage, 2),
        "with_stats": with_stats,
        "outliers": outliers,
    }

    return {
        "score": round(score, 2),
        "status": status,
        "issues": issues,
        "metrics": metrics,
    }

