from __future__ import annotations

from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import linregress
from sklearn.cluster import KMeans


def detect_trend(values: Iterable[float], timepoints: Iterable[float]) -> Tuple[float, float, str]:
    vals = np.array(list(values), dtype=float)
    times = np.array(list(timepoints), dtype=float)
    if len(vals) != len(times) or len(vals) < 2:
        return np.nan, np.nan, "stable"

    mask = ~(np.isnan(vals) | np.isnan(times))
    vals = vals[mask]
    times = times[mask]
    if len(vals) < 2:
        return np.nan, np.nan, "stable"

    slope, _, rvalue, pvalue, _ = linregress(times, vals)
    direction = "stable"
    if not np.isnan(pvalue) and pvalue < 0.05:
        if slope > 0:
            direction = "increasing"
        elif slope < 0:
            direction = "decreasing"
    return slope, pvalue, direction


def cluster_trajectories(df: pd.DataFrame, timepoints: List[float], n_clusters: int = 3) -> np.ndarray:
    if df.empty or df.shape[0] < n_clusters:
        return np.zeros(df.shape[0], dtype=int)
    data = df.to_numpy(dtype=float)
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
    labels = kmeans.fit_predict(data)
    return labels


def find_temporal_patterns(df: pd.DataFrame, timepoints: List[float], n_clusters: int = 3) -> pd.DataFrame:
    """Analyze temporal trends and cluster feature trajectories.

    Args:
        df: DataFrame with rows=features, columns=timepoints (numeric values)
        timepoints: list of time values aligned to df columns
        n_clusters: number of clusters for trajectory grouping

    Returns:
        DataFrame with feature, slope, pvalue, trend_direction, cluster_id
    """
    features = df.index.tolist()
    slopes: List[float] = []
    pvals: List[float] = []
    directions: List[str] = []

    for _, row in df.iterrows():
        s, p, d = detect_trend(row.values, timepoints)
        slopes.append(s)
        pvals.append(p)
        directions.append(d)

    cluster_ids = cluster_trajectories(df, timepoints, n_clusters=n_clusters)

    return pd.DataFrame(
        {
            "feature": features,
            "slope": slopes,
            "pvalue": pvals,
            "trend_direction": directions,
            "cluster_id": cluster_ids,
        }
    )


__all__ = ["detect_trend", "cluster_trajectories", "find_temporal_patterns"]
