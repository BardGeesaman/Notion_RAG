from __future__ import annotations

import numpy as np
import pytest


pytest.importorskip("sklearn")


def _make_data(seed: int = 0, n_samples: int = 120, n_features: int = 10):
    rng = np.random.RandomState(seed)
    X = rng.normal(size=(n_samples, n_features))
    score = 1.5 * X[:, 0] - 1.0 * X[:, 1] + 0.25 * rng.normal(size=(n_samples,))
    y = (score > np.median(score)).astype(int)
    names = [f"f{i}" for i in range(n_features)]
    return X, y, names


def test_cv_feature_importance_fit_returns_rankings(monkeypatch):
    from amprenta_rag.ml.biomarker import importance as imod

    # Encourage determinism if sklearn is parallelized: force n_jobs=1.
    from sklearn.ensemble import RandomForestClassifier as _RFC  # type: ignore

    def rfc_factory(*, n_estimators, random_state, n_jobs=-1, **kwargs):  # noqa: ANN001
        _ = n_jobs
        return _RFC(n_estimators=n_estimators, random_state=random_state, n_jobs=1, **kwargs)

    monkeypatch.setattr(imod, "RandomForestClassifier", rfc_factory)

    from amprenta_rag.ml.biomarker.importance import CVFeatureImportance

    X, y, names = _make_data(seed=0)
    imp = CVFeatureImportance(n_folds=4, n_estimators=50)
    imp.fit(X, y, feature_names=names, random_state=0)
    ranked = imp.get_ranked_features()

    assert len(ranked) == len(names)
    means = [m for _, m, _ in ranked]
    assert all(means[i] >= means[i + 1] for i in range(len(means) - 1))


def test_cv_feature_importance_reproducible_with_seed(monkeypatch):
    from amprenta_rag.ml.biomarker import importance as imod

    from sklearn.ensemble import RandomForestClassifier as _RFC  # type: ignore

    def rfc_factory(*, n_estimators, random_state, n_jobs=-1, **kwargs):  # noqa: ANN001
        _ = n_jobs
        return _RFC(n_estimators=n_estimators, random_state=random_state, n_jobs=1, **kwargs)

    monkeypatch.setattr(imod, "RandomForestClassifier", rfc_factory)

    from amprenta_rag.ml.biomarker.importance import CVFeatureImportance

    X, y, names = _make_data(seed=3)

    imp1 = CVFeatureImportance(n_folds=4, n_estimators=50)
    imp1.fit(X, y, feature_names=names, random_state=123)
    r1 = imp1.get_ranked_features()

    imp2 = CVFeatureImportance(n_folds=4, n_estimators=50)
    imp2.fit(X, y, feature_names=names, random_state=123)
    r2 = imp2.get_ranked_features()

    top1 = [x[0] for x in r1[:5]]
    top2 = [x[0] for x in r2[:5]]
    assert top1 == top2
    assert np.allclose([x[1] for x in r1], [x[1] for x in r2], atol=1e-6)


