from __future__ import annotations

import numpy as np
import pytest

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


def test_train_target_model_returns_metrics(monkeypatch):
    from amprenta_rag.ml.qsar import trainer as tmod

    # Mock dataset loader output (enough rows for splits + calibrator >=10).
    X = np.zeros((120, 2054), dtype=np.float32)
    y = np.array([1] * 40 + [0] * 80, dtype=np.int64)
    smiles = ["C"] * 120
    meta = {"source": "chembl"}

    class FakeLoader:
        def load_target(self, *args, **kwargs):  # noqa: ANN001
            return X, y, smiles, meta

    class FakeEnsemble:
        def __init__(self, n_models=5, base_params=None):  # noqa: ANN001
            self.n_models = n_models
            self.base_params = dict(base_params or {})
            self._trained = False
            self.training_centroid = None
            self.feature_dim = None
            self.task = "classification"

        def fit(self, Xtr, ytr):  # noqa: ANN001
            self._trained = True
            self.feature_dim = int(Xtr.shape[1])
            self.training_centroid = np.asarray(Xtr).mean(axis=0).astype(np.float32)
            return self

        def predict_proba(self, Xte):  # noqa: ANN001
            n = int(np.asarray(Xte).shape[0])
            mean = np.linspace(0.1, 0.9, n)
            std = np.zeros((n,), dtype=float)
            return mean, std

        def to_artifact(self):  # noqa: ANN001
            return {
                "models": [],
                "training_centroid": self.training_centroid,
                "n_models": int(self.n_models),
                "feature_dim": int(self.feature_dim or 2054),
                "task": "classification",
            }

    class FakeCalibrator:
        def __init__(self, method="isotonic"):  # noqa: ANN001
            self.method = method
            self.calibrator = object()

        def fit(self, y_pred, y_true):  # noqa: ANN001
            return self

        def calibrate(self, y_pred):  # noqa: ANN001
            return np.asarray(y_pred, dtype=float)

        def compute_ece(self, y_pred, y_true, n_bins=10):  # noqa: ANN001
            return 0.0

    monkeypatch.setattr(tmod, "TargetDatasetLoader", lambda: FakeLoader())
    monkeypatch.setattr(tmod, "BootstrapEnsemble", FakeEnsemble)
    monkeypatch.setattr(tmod, "CalibrationWrapper", FakeCalibrator)

    out = tmod.train_target_model(target="EGFR", register=False)
    assert "metrics" in out
    assert set(out["metrics"].keys()) >= {"auc", "accuracy", "ece"}
    assert out["model_name"] == "qsar_EGFR_ensemble"


def test_train_target_model_registers_model(monkeypatch):
    from amprenta_rag.ml.qsar import trainer as tmod

    X = np.zeros((120, 2054), dtype=np.float32)
    y = np.array([1] * 40 + [0] * 80, dtype=np.int64)

    class FakeLoader:
        def load_target(self, *args, **kwargs):  # noqa: ANN001
            return X, y, ["C"] * 120, {}

    class FakeEnsemble:
        def __init__(self, n_models=5, base_params=None):  # noqa: ANN001
            self.base_params = dict(base_params or {})

        def fit(self, Xtr, ytr):  # noqa: ANN001
            return self

        def predict_proba(self, Xte):  # noqa: ANN001
            n = int(np.asarray(Xte).shape[0])
            return np.full((n,), 0.5, dtype=float), np.zeros((n,), dtype=float)

        def to_artifact(self):  # noqa: ANN001
            return {"models": [], "training_centroid": np.zeros((2054,), dtype=np.float32), "n_models": 5, "feature_dim": 2054, "task": "classification"}

    class FakeCalibrator:
        def __init__(self, method="isotonic"):  # noqa: ANN001
            self.calibrator = object()

        def fit(self, y_pred, y_true):  # noqa: ANN001
            return self

        def calibrate(self, y_pred):  # noqa: ANN001
            return np.asarray(y_pred, dtype=float)

        def compute_ece(self, y_pred, y_true, n_bins=10):  # noqa: ANN001
            return 0.0

    class FakeRegistry:
        def __init__(self):
            self.calls = []

        def register_model(self, **kwargs):  # noqa: ANN001
            self.calls.append(kwargs)
            return object()

    reg = FakeRegistry()

    monkeypatch.setattr(tmod, "TargetDatasetLoader", lambda: FakeLoader())
    monkeypatch.setattr(tmod, "BootstrapEnsemble", FakeEnsemble)
    monkeypatch.setattr(tmod, "CalibrationWrapper", FakeCalibrator)
    monkeypatch.setattr(tmod, "get_registry", lambda: reg)

    out = tmod.train_target_model(target="EGFR", register=True)
    assert reg.calls, "Expected register_model to be called"
    call = reg.calls[0]
    assert call["name"] == "qsar_EGFR_ensemble"
    assert call["model_type"] == "qsar_classification"
    assert out["model_name"] == "qsar_EGFR_ensemble"


def test_class_imbalance_handled(monkeypatch):
    from amprenta_rag.ml.qsar import trainer as tmod

    X = np.zeros((120, 2054), dtype=np.float32)
    # 10% positives, 90% negatives
    y = np.array([1] * 12 + [0] * 108, dtype=np.int64)

    class FakeLoader:
        def load_target(self, *args, **kwargs):  # noqa: ANN001
            return X, y, ["C"] * 120, {}

    captured = {}

    class FakeEnsemble:
        def __init__(self, n_models=5, base_params=None):  # noqa: ANN001
            captured["base_params"] = dict(base_params or {})

        def fit(self, Xtr, ytr):  # noqa: ANN001
            return self

        def predict_proba(self, Xte):  # noqa: ANN001
            n = int(np.asarray(Xte).shape[0])
            return np.full((n,), 0.5, dtype=float), np.zeros((n,), dtype=float)

        def to_artifact(self):  # noqa: ANN001
            return {"models": [], "training_centroid": np.zeros((2054,), dtype=np.float32), "n_models": 5, "feature_dim": 2054, "task": "classification"}

    class FakeCalibrator:
        def __init__(self, method="isotonic"):  # noqa: ANN001
            self.calibrator = object()

        def fit(self, y_pred, y_true):  # noqa: ANN001
            return self

        def calibrate(self, y_pred):  # noqa: ANN001
            return np.asarray(y_pred, dtype=float)

        def compute_ece(self, y_pred, y_true, n_bins=10):  # noqa: ANN001
            return 0.0

    monkeypatch.setattr(tmod, "TargetDatasetLoader", lambda: FakeLoader())
    monkeypatch.setattr(tmod, "BootstrapEnsemble", FakeEnsemble)
    monkeypatch.setattr(tmod, "CalibrationWrapper", FakeCalibrator)

    _ = tmod.train_target_model(target="EGFR", register=False)
    assert "scale_pos_weight" in captured["base_params"]
    # scale_pos_weight computed on train split; should be > 1 for imbalanced data.
    assert float(captured["base_params"]["scale_pos_weight"]) > 1.0


