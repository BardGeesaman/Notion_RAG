from __future__ import annotations

import numpy as np
import pytest


pytest.importorskip("sklearn")


def _make_data(seed: int = 0, n_samples: int = 80, n_features: int = 12):
    rng = np.random.RandomState(seed)
    X = rng.normal(size=(n_samples, n_features))
    # Make y correlated with feature 0
    score = 2.0 * X[:, 0] + 0.25 * rng.normal(size=(n_samples,))
    y = (score > np.median(score)).astype(int)
    feature_names = [f"f{i}" for i in range(n_features)]
    return X, y, feature_names


def test_stability_selector_fit_returns_rankings():
    from amprenta_rag.ml.biomarker.stability import StabilitySelector

    X, y, names = _make_data(seed=0)
    sel = StabilitySelector(n_bootstrap=25, threshold=0.0)
    sel.fit(X, y, feature_names=names, random_state=0)
    ranked = sel.get_ranked_features()

    assert isinstance(ranked, list)
    assert len(ranked) == len(names)
    assert ranked[0][1] >= ranked[-1][1]
    # Sorted descending by selection frequency
    freqs = [r[1] for r in ranked]
    assert all(freqs[i] >= freqs[i + 1] for i in range(len(freqs) - 1))


def test_stability_selector_threshold_filtering():
    from amprenta_rag.ml.biomarker.stability import StabilitySelector

    X, y, names = _make_data(seed=1)

    low = StabilitySelector(n_bootstrap=25, threshold=0.0)
    low.fit(X, y, feature_names=names, random_state=0)
    low_ranked = low.get_ranked_features()

    high = StabilitySelector(n_bootstrap=25, threshold=0.8)
    high.fit(X, y, feature_names=names, random_state=0)
    high_ranked = high.get_ranked_features()

    assert len(high_ranked) <= len(low_ranked)
    assert all(freq >= 0.8 for _, freq, _ in high_ranked)


