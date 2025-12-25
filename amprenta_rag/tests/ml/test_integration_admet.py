from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
from fastapi.testclient import TestClient

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


class _DummyEnsembleClassifier:
    def __init__(self, p: float):
        self.p = float(p)

    def predict_proba(self, X):  # noqa: ANN001
        n = len(X)
        p1 = np.full((n,), self.p, dtype=float)
        return np.column_stack([1.0 - p1, p1])


def test_full_prediction_pipeline(monkeypatch):
    from amprenta_rag.ml.admet.predictor import ADMETPredictor

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
        "metrics": {"ece": 0.2},
        "trained_at": "2025-12-24T00:00:00Z",
    }

    monkeypatch.setattr(predictor.registry, "get_active_model", lambda name: SimpleNamespace(id="x", version="1.0.0"))
    monkeypatch.setattr(predictor.registry, "load_model", lambda model_id: artifact)  # noqa: ARG005

    out = predictor.predict_with_uncertainty(["CCO"], endpoints=["herg"])
    assert isinstance(out, list) and out
    r0 = out[0]
    assert "predictions" in r0
    p0 = r0["predictions"]["herg"]
    assert p0["ci_low"] < p0["mean"] < p0["ci_high"]


def test_api_endpoint_returns_valid_response(monkeypatch):
    from fastapi import FastAPI
    from amprenta_rag.api.schemas import ADMETPredictResponse
    import amprenta_rag.api.routers.admet as admet_router

    # Build a minimal app to avoid importing the full API (which may require optional deps like python-multipart).
    app = FastAPI()
    app.include_router(admet_router.router, prefix="/api")

    # Override DB dependency to avoid Postgres.
    from amprenta_rag.database.base import get_db

    def _dummy_get_db():
        yield None

    app.dependency_overrides[get_db] = _dummy_get_db

    class _Pred:
        def predict_with_uncertainty(self, smiles_list, endpoints=None):  # noqa: ANN001
            return [
                {
                    "smiles": smiles_list[0],
                    "predictions": {
                        "herg": {
                            "mean": 0.5,
                            "std": 0.1,
                            "ci_low": 0.3,
                            "ci_high": 0.7,
                            "in_domain": True,
                            "similarity": 0.9,
                            "calibrated": False,
                        }
                    },
                    "error": None,
                }
            ]

        def predict(self, smiles_list, endpoints=None, include_shap=False):  # noqa: ANN001
            return [{"smiles": smiles_list[0], "herg": {"probability": 0.5, "class": 0}}]

    monkeypatch.setattr(admet_router, "get_admet_predictor", lambda: _Pred())
    monkeypatch.setattr(admet_router, "get_registry", lambda: SimpleNamespace(get_active_model=lambda name: None))

    client = TestClient(app)
    resp = client.post("/api/admet/predict", json={"smiles": ["CCO"], "endpoints": ["herg"], "include_uncertainty": True})
    assert resp.status_code == 200
    body = resp.json()
    ADMETPredictResponse.model_validate(body)


def test_calibration_improves_ece():
    from amprenta_rag.ml.admet.calibration import CalibrationWrapper, reliability_diagram

    y_pred = np.linspace(0.0, 1.0, 200)
    # Strongly miscalibrated: make outcomes too "peaky"
    y_true = (y_pred > 0.8).astype(float)
    ece_before = float(reliability_diagram(y_pred, y_true, n_bins=10)["ece"])

    cal = CalibrationWrapper(method="isotonic").fit(y_pred, y_true)
    y_cal = cal.calibrate(y_pred)
    ece_after = float(reliability_diagram(y_cal, y_true, n_bins=10)["ece"])

    assert ece_after <= ece_before + 1e-9
    assert ece_after < ece_before


def test_applicability_widens_ci_for_ood():
    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker

    chk = ApplicabilityChecker(threshold=0.3)
    chk.training_centroid = np.ones((10,), dtype=np.float32)

    # identical features -> in-domain
    X_in = np.ones((1, 10), dtype=np.float32)
    in_domain_in, _ = chk.check(X_in)

    # dissimilar -> OOD
    X_out = np.zeros((1, 10), dtype=np.float32)
    in_domain_out, _ = chk.check(X_out)

    ci_low = np.array([0.0], dtype=float)
    ci_high = np.array([1.0], dtype=float)

    low_in, high_in = chk.widen_ci(ci_low, ci_high, in_domain_in)
    low_out, high_out = chk.widen_ci(ci_low, ci_high, in_domain_out)

    assert float(high_in - low_in) == 1.0
    assert float(high_out - low_out) == 2.0


