from __future__ import annotations

import numpy as np
import pandas as pd

from amprenta_rag.analysis import timeseries


def test_detect_trend_increasing():
    slope, p, direction = timeseries.detect_trend([1, 2, 3, 4], [0, 1, 2, 3])
    assert direction in {"increasing", "stable"}
    assert not np.isnan(slope)


def test_detect_trend_insufficient_or_nans():
    slope, p, direction = timeseries.detect_trend([np.nan, 1], [0, 1])
    assert np.isnan(slope)
    assert direction == "stable"


def test_cluster_trajectories_handles_small():
    df = pd.DataFrame()
    labels = timeseries.cluster_trajectories(df, [], n_clusters=3)
    assert labels.size == 0


def test_find_temporal_patterns_outputs_dataframe():
    data = pd.DataFrame([[1, 2, 3], [3, 2, 1]], index=["f1", "f2"])
    tp = [0, 1, 2]
    out = timeseries.find_temporal_patterns(data, tp, n_clusters=2)
    assert set(out.columns) == {"feature", "slope", "pvalue", "trend_direction", "cluster_id"}
    assert len(out) == 2

