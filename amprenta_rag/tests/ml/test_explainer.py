from __future__ import annotations

import numpy as np
import pytest


shap = pytest.importorskip("shap")
xgb = pytest.importorskip("xgboost")
_ = (shap, xgb)


from amprenta_rag.ml.admet.explainer import EnsembleSHAPExplainer  # noqa: E402


class _DummyTreeExplainer:
    def __init__(self, model):
        self.model = model
        self.expected_value = float(getattr(model, "base", 0.0))

    def shap_values(self, X):
        # Deterministic SHAP pattern based on model.idx
        X = np.asarray(X)
        sv = np.zeros((X.shape[0], 2054), dtype=float)
        i = int(getattr(self.model, "idx", 0))
        sv[:, i] = 10.0 + i
        sv[:, (i + 1) % 2054] = -5.0
        return sv


class _M:
    def __init__(self, idx: int, base: float = 1.0):
        self.idx = idx
        self.base = base


def _artifact():
    models = [_M(i, base=1.0 + i) for i in range(5)]
    return {"ensemble": {"models": models}}


def test_explain_prediction_returns_top_features(monkeypatch):
    ex = EnsembleSHAPExplainer(_artifact())
    monkeypatch.setattr(ex._shap, "TreeExplainer", _DummyTreeExplainer)

    X = np.zeros((1, 2054), dtype=float)
    out = ex.explain_prediction(X, top_k=3)

    assert out["shap_values"].shape == (2054,)
    assert len(out["top_features"]) == 3
    assert [t["rank"] for t in out["top_features"]] == [1, 2, 3]


def test_global_importance_sorted_descending(monkeypatch):
    ex = EnsembleSHAPExplainer(_artifact())
    monkeypatch.setattr(ex._shap, "TreeExplainer", _DummyTreeExplainer)

    X = np.zeros((4, 2054), dtype=float)
    out = ex.compute_global_importance(X)

    assert len(out) == 2054
    imps = [d["importance"] for d in out[:20]]
    assert all(imps[i] >= imps[i + 1] for i in range(len(imps) - 1))


def test_other_sum_calculation(monkeypatch):
    ex = EnsembleSHAPExplainer(_artifact())
    monkeypatch.setattr(ex._shap, "TreeExplainer", _DummyTreeExplainer)

    X = np.zeros((1, 2054), dtype=float)
    out = ex.explain_prediction(X, top_k=2)

    sv = out["shap_values"]
    top_idx = {int(t["name"].split("_")[-1]) if t["name"].startswith("MorganBit_") else None for t in out["top_features"]}
    # 'other_sum' should match sum of sv excluding the top two features.
    # (We can't reliably parse descriptor names here; for this synthetic test, top features are Morgan bits.)
    top_idx = {i for i in top_idx if i is not None}
    expected = float(np.sum([sv[i] for i in range(2054) if i not in top_idx]))
    assert abs(out["other_sum"] - expected) < 1e-8


