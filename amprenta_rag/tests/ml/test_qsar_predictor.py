from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict
from uuid import UUID, uuid4

import numpy as np
import pytest

pytest.importorskip("rdkit")


@dataclass
class FakeMLModel:
    id: UUID
    name: str
    version: str = "1.0.0"
    metrics: Dict[str, Any] = None  # type: ignore[assignment]


class FakeXGB:
    def __init__(self, p: float):
        self.p = float(p)

    def predict_proba(self, X):  # noqa: ANN001
        n = int(np.asarray(X).shape[0])
        p1 = np.full((n,), self.p, dtype=float)
        p0 = 1.0 - p1
        return np.vstack([p0, p1]).T


class FakeCalibrator:
    def calibrate(self, y_pred):  # noqa: ANN001
        return np.clip(np.asarray(y_pred, dtype=float), 0.0, 1.0)


def test_list_available_targets(monkeypatch):
    from amprenta_rag.ml.qsar.predictor import TargetQSARPredictor

    class FakeRegistry:
        def list_models(self):  # noqa: ANN001
            return [
                FakeMLModel(id=uuid4(), name="qsar_EGFR_ensemble", metrics={"auc": 0.8}),
                FakeMLModel(id=uuid4(), name="admet_herg_ensemble", metrics={}),
                FakeMLModel(id=uuid4(), name="qsar_BRAF_ensemble", metrics={"auc": 0.7}),
            ]

    pred = TargetQSARPredictor()
    pred.registry = FakeRegistry()  # type: ignore[assignment]
    out = pred.list_available_targets()
    targets = {d["target"] for d in out}
    assert targets == {"EGFR", "BRAF"}


def test_predict_returns_probabilities(monkeypatch):
    from amprenta_rag.ml.qsar.predictor import TargetQSARPredictor
    from amprenta_rag.ml.admet import predictor as predmod

    # Feature extraction should not require RDKit in this unit test.
    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", lambda self, s: np.zeros((2054,), dtype=np.float32))

    egfr_model = FakeMLModel(id=uuid4(), name="qsar_EGFR_ensemble", metrics={"auc": 0.9})

    artifact = {
        "ensemble": {"models": [FakeXGB(0.8), FakeXGB(0.6)], "task": "classification"},
        "calibrator": FakeCalibrator(),
        "applicability": {"threshold": 0.3, "centroid": np.zeros((2054,), dtype=np.float32)},
        "metadata": {"target": "EGFR", "compound_count": 123},
    }

    class FakeRegistry:
        def get_active_model(self, name: str):  # noqa: ANN001
            return egfr_model if name == "qsar_EGFR_ensemble" else None

        def load_model(self, model_id: UUID):  # noqa: ANN001
            assert model_id == egfr_model.id
            return artifact

    pred = TargetQSARPredictor()
    pred.registry = FakeRegistry()  # type: ignore[assignment]

    out = pred.predict(["CCO"], targets=["EGFR"])
    assert len(out) == 1
    p = out[0]["predictions"]["EGFR"]
    assert 0.0 <= p["probability"] <= 1.0
    assert p["ci_low"] <= p["probability"] <= p["ci_high"]
    assert p["calibrated"] is True
    assert isinstance(p["in_domain"], bool)


def test_predict_invalid_smiles(monkeypatch):
    from amprenta_rag.ml.qsar.predictor import TargetQSARPredictor
    from amprenta_rag.ml.admet import predictor as predmod

    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", lambda self, s: None)

    class FakeRegistry:
        def get_active_model(self, name: str):  # noqa: ANN001
            return FakeMLModel(id=uuid4(), name=name)

        def load_model(self, model_id: UUID):  # noqa: ANN001
            return {"ensemble": {"models": [FakeXGB(0.5)], "task": "classification"}, "applicability": {"threshold": 0.3, "centroid": np.zeros((2054,), dtype=np.float32)}}

    pred = TargetQSARPredictor()
    pred.registry = FakeRegistry()  # type: ignore[assignment]
    out = pred.predict(["NOT_A_SMILES"], targets=["EGFR"])
    assert out[0]["error"] == "Invalid SMILES"


def test_predict_missing_model_handled(monkeypatch):
    from amprenta_rag.ml.qsar.predictor import TargetQSARPredictor
    from amprenta_rag.ml.admet import predictor as predmod

    monkeypatch.setattr(predmod.ADMETPredictor, "_get_features", lambda self, s: np.zeros((2054,), dtype=np.float32))

    class FakeRegistry:
        def get_active_model(self, name: str):  # noqa: ANN001
            return None

    pred = TargetQSARPredictor()
    pred.registry = FakeRegistry()  # type: ignore[assignment]
    out = pred.predict(["CCO"], targets=["EGFR"])
    assert out[0]["predictions"] == {}
    assert out[0]["error"] == "No target models available"


