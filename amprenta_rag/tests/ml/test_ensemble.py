from __future__ import annotations

import numpy as np
import pytest

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


class _DummyXGBClassifier:
    def __init__(self, **kwargs):  # noqa: ANN003
        self.seed = int(kwargs.get("random_state", 0))
        self.bias = 0.0

    def fit(self, X, y):  # noqa: ANN001
        # Fit a simple bias so different bootstrap samples yield different outputs.
        y = np.asarray(y, dtype=float)
        self.bias = float(np.mean(y)) if len(y) else 0.0
        return self

    def predict_proba(self, X):  # noqa: ANN001
        X = np.asarray(X, dtype=float)
        rng = np.random.default_rng(self.seed)
        w = rng.normal(size=(X.shape[1],))
        logits = X @ w / max(1.0, float(X.shape[1])) + self.bias
        p1 = 1.0 / (1.0 + np.exp(-logits))
        p1 = np.clip(p1, 0.0, 1.0)
        return np.column_stack([1.0 - p1, p1])


class _DummyXGBRegressor:
    def __init__(self, **kwargs):  # noqa: ANN003
        self.seed = int(kwargs.get("random_state", 0))
        self.bias = 0.0

    def fit(self, X, y):  # noqa: ANN001
        y = np.asarray(y, dtype=float)
        self.bias = float(np.mean(y)) if len(y) else 0.0
        return self

    def predict(self, X):  # noqa: ANN001
        X = np.asarray(X, dtype=float)
        rng = np.random.default_rng(self.seed)
        w = rng.normal(size=(X.shape[1],))
        return (X @ w / max(1.0, float(X.shape[1])) + self.bias).astype(float)


def test_bootstrap_sampling_diversity(monkeypatch):
    from amprenta_rag.ml.admet import ensemble as emod

    monkeypatch.setattr(emod, "_require_xgboost", lambda: (_DummyXGBClassifier, _DummyXGBRegressor))

    X = np.random.default_rng(0).integers(0, 2, size=(100, 32)).astype(np.float32)
    y = (X[:, 0] > 0).astype(int)

    ens = emod.BootstrapEnsemble(n_models=5, base_params={"n_estimators": 5}).fit(X, y)
    idxs = [tuple(i.tolist()) for i in ens.bootstrap_indices_]
    assert len(set(idxs)) > 1, "Expected bootstrap samples to differ across models"


def test_ensemble_predict_returns_mean_std(monkeypatch):
    from amprenta_rag.ml.admet import ensemble as emod

    monkeypatch.setattr(emod, "_require_xgboost", lambda: (_DummyXGBClassifier, _DummyXGBRegressor))

    X = np.random.default_rng(1).normal(size=(120, 16)).astype(np.float32)
    y = np.random.default_rng(2).normal(size=(120,)).astype(np.float32)

    ens = emod.BootstrapEnsemble(n_models=3, base_params={"n_estimators": 5}).fit(X, y)
    mean, std = ens.predict(X[:10])
    assert mean.shape == (10,)
    assert std.shape == (10,)
    assert np.all(std >= 0.0)


def test_ensemble_predict_proba_classification(monkeypatch):
    from amprenta_rag.ml.admet import ensemble as emod

    monkeypatch.setattr(emod, "_require_xgboost", lambda: (_DummyXGBClassifier, _DummyXGBRegressor))

    X = np.random.default_rng(3).integers(0, 2, size=(200, 24)).astype(np.float32)
    y = (X[:, 1] + X[:, 2] > 0).astype(int)

    ens = emod.BootstrapEnsemble(n_models=4, base_params={"n_estimators": 5}).fit(X, y)
    mean_p, std_p = ens.predict_proba(X[:15])
    assert mean_p.shape == (15,)
    assert std_p.shape == (15,)
    assert np.all((mean_p >= 0.0) & (mean_p <= 1.0))
    assert np.all(std_p >= 0.0)


def test_pairwise_correlation_warning(caplog):
    from amprenta_rag.ml.admet.ensemble import BootstrapEnsemble

    preds = np.array(
        [
            [0.1, 0.2, 0.3, 0.4],
            [0.11, 0.21, 0.31, 0.41],
            [0.09, 0.19, 0.29, 0.39],
        ],
        dtype=float,
    )

    with caplog.at_level("WARNING"):
        BootstrapEnsemble._warn_if_high_corr(preds, threshold=0.95)

    assert any("High pairwise ensemble correlation" in r.message for r in caplog.records)


