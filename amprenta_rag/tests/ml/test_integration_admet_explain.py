from __future__ import annotations

from types import SimpleNamespace

import numpy as np


class _DummyEnsembleClassifier:
    def __init__(self, p: float):
        self.p = float(p)

    def predict_proba(self, X):  # noqa: ANN001
        n = len(X)
        p1 = np.full((n,), self.p, dtype=float)
        return np.column_stack([1.0 - p1, p1])


def test_full_explain_pipeline(monkeypatch):
    """Mock registry + mock explainer and verify predictor returns shap payload when include_shap=True."""
    from amprenta_rag.ml.admet.predictor import ADMETPredictor
    import amprenta_rag.ml.admet.explainer as expmod

    predictor = ADMETPredictor()
    monkeypatch.setattr(predictor, "_get_features", lambda smiles: np.ones((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.2), _DummyEnsembleClassifier(0.4)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": None,
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
    }

    monkeypatch.setattr(
        predictor.registry,
        "get_active_model",
        lambda name: SimpleNamespace(id="x", version="1.0.0"),
    )
    monkeypatch.setattr(predictor.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    class _DummyExplainer:
        def __init__(self, _artifact):  # noqa: ANN001
            pass

        def explain_prediction(self, X, top_k: int = 10):  # noqa: ANN001
            assert np.asarray(X).shape == (1, 2054)
            assert top_k == 7
            return {
                "shap_values": np.zeros((2054,), dtype=float),
                "top_features": [{"name": "MorganBit_0000", "value": 1.0, "rank": 1}],
                "base_value": 0.5,
                "other_sum": 0.0,
            }

    monkeypatch.setattr(expmod, "EnsembleSHAPExplainer", _DummyExplainer)

    out = predictor.predict_with_uncertainty(["CCO"], endpoints=["herg"], include_shap=True, shap_top_k=7)
    pred = out[0]["predictions"]["herg"]
    assert pred["mean"] is not None
    assert pred["shap"]["top_features"][0]["name"] == "MorganBit_0000"
    assert pred["shap"]["base_value"] == 0.5


def test_global_importance_accumulation():
    """MVP global importance: mean abs(value) over accumulated top_features from multiple explains."""

    hist = [
        {"top_features": [{"name": "MorganBit_0001", "value": 2.0, "rank": 1}]},
        {"top_features": [{"name": "MorganBit_0001", "value": -4.0, "rank": 1}]},
        {"top_features": [{"name": "MolWt", "value": 1.0, "rank": 1}]},
    ]

    acc: dict[str, float] = {}
    n = 0
    for item in hist:
        top = item.get("top_features")
        if not isinstance(top, list):
            continue
        n += 1
        for f in top:
            name = str(f.get("name") or "")
            val = float(f.get("value", 0.0))
            acc[name] = acc.get(name, 0.0) + abs(val)

    imp = [{"name": k, "importance": v / float(n)} for k, v in acc.items()]
    imp.sort(key=lambda x: float(x["importance"]), reverse=True)

    # MorganBit_0001: (|2|+|4|)/3 = 2.0
    assert imp[0]["name"] == "MorganBit_0001"
    assert abs(float(imp[0]["importance"]) - 2.0) < 1e-9


