from __future__ import annotations

import numpy as np


class _DummyEnsembleClassifier:
    def __init__(self, p: float):
        self.p = float(p)

    def predict_proba(self, X):  # noqa: ANN001
        n = len(X)
        p1 = np.full((n,), self.p, dtype=float)
        return np.column_stack([1.0 - p1, p1])


class _DummyEnsembleRegressor:
    def __init__(self, v: float):
        self.v = float(v)

    def predict(self, X):  # noqa: ANN001
        return np.full((len(X),), self.v, dtype=float)


class _DummyCalibrator:
    def calibrate(self, y_pred):  # noqa: ANN001
        y_pred = np.asarray(y_pred, dtype=float)
        return np.clip(y_pred * 0.5, 0.0, 1.0)


def test_predict_with_uncertainty_returns_ci(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.ones((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.2), _DummyEnsembleClassifier(0.4)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": _DummyCalibrator(),
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
        "trained_at": "2025-12-24T00:00:00Z",
        "metrics": {"ece": 0.1},
    }

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    out = p.predict_with_uncertainty(["CCO"], endpoints=["herg"])
    pred = out[0]["predictions"]["herg"]
    assert pred["ci_low"] < pred["mean"] < pred["ci_high"]
    assert pred["calibrated"] is True


def test_predict_with_uncertainty_ood_warning(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.zeros((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.5), _DummyEnsembleClassifier(0.5)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": None,
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
    }

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    out = p.predict_with_uncertainty(["CCO"], endpoints=["herg"])
    pred = out[0]["predictions"]["herg"]
    assert pred["in_domain"] is False
    # OOD widening should expand CI beyond the base +/- 1.96*std (std=0 here so CI stays same)
    assert pred["similarity"] is not None
    assert pred["similarity"] <= 1e-6


def test_predict_with_uncertainty_backward_compat(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.zeros((2054,), dtype=np.float32))

    class _DummyClassifier:
        def predict_proba(self, X):  # noqa: ANN001
            return np.tile(np.array([[0.2, 0.8]], dtype=float), (len(X), 1))

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: _DummyClassifier())  # noqa: ARG005

    out = p.predict(["CCO"], endpoints=["herg"], include_shap=False)
    assert out[0]["herg"]["probability"] == 0.8


def test_predict_with_uncertainty_invalid_smiles(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: None)

    out = p.predict_with_uncertainty(["NOT_A_SMILES"], endpoints=["herg"])
    assert out[0]["error"] is not None


def test_predict_with_uncertainty_shap_disabled(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.ones((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.2), _DummyEnsembleClassifier(0.4)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": _DummyCalibrator(),
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
    }

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    out = p.predict_with_uncertainty(["CCO"], endpoints=["herg"])
    pred = out[0]["predictions"]["herg"]
    assert "shap" not in pred


def test_predict_with_uncertainty_shap_included(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod
    import amprenta_rag.ml.admet.explainer as expmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.ones((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.2), _DummyEnsembleClassifier(0.4)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": _DummyCalibrator(),
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
    }

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    called: dict[str, int] = {"k": 0}

    class _DummyExplainer:
        def __init__(self, _artifact):  # noqa: ANN001
            pass

        def explain_prediction(self, X, top_k: int = 10):  # noqa: ANN001
            assert np.asarray(X).shape == (1, 2054)
            called["k"] = int(top_k)
            return {
                "shap_values": np.zeros((2054,), dtype=float),
                "top_features": [{"name": "MorganBit_0000", "value": 1.0, "rank": 1}],
                "base_value": 0.5,
                "other_sum": 0.0,
            }

    monkeypatch.setattr(expmod, "EnsembleSHAPExplainer", _DummyExplainer)

    out = p.predict_with_uncertainty(["CCO"], endpoints=["herg"], include_shap=True, shap_top_k=7)
    pred = out[0]["predictions"]["herg"]
    assert "shap" in pred
    assert pred["shap"]["top_features"][0]["name"] == "MorganBit_0000"
    assert pred["shap"]["base_value"] == 0.5
    assert pred["shap"]["other_sum"] == 0.0
    assert called["k"] == 7


def test_predict_with_uncertainty_shap_graceful_degradation(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod
    import amprenta_rag.ml.admet.explainer as expmod

    p = predmod.ADMETPredictor()
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.ones((2054,), dtype=np.float32))

    artifact = {
        "ensemble": {
            "models": [_DummyEnsembleClassifier(0.2), _DummyEnsembleClassifier(0.4)],
            "training_centroid": np.ones((2054,), dtype=np.float32),
            "n_models": 2,
            "feature_dim": 2054,
            "task": "classification",
        },
        "calibrator": _DummyCalibrator(),
        "applicability": {"centroid": np.ones((2054,), dtype=np.float32), "threshold": 0.3},
    }

    class _Rec:
        def __init__(self):
            self.id = "dummy"

    monkeypatch.setattr(p.registry, "get_active_model", lambda name: _Rec())
    monkeypatch.setattr(p.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    def _boom(_artifact):  # noqa: ANN001
        raise ImportError("no shap")

    monkeypatch.setattr(expmod, "EnsembleSHAPExplainer", _boom)

    out = p.predict_with_uncertainty(["CCO"], endpoints=["herg"], include_shap=True)
    pred = out[0]["predictions"]["herg"]
    assert pred["shap"]["error"] == "shap not installed"


