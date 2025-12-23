from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4

import numpy as np
import pytest
from fastapi.testclient import TestClient


class _DummyClassifier:
    def predict_proba(self, X):  # noqa: ANN001
        # return class 1 prob = 0.8
        return np.tile(np.array([[0.2, 0.8]], dtype=float), (len(X), 1))


class _DummyRegressor:
    def predict(self, X):  # noqa: ANN001
        return np.full((len(X),), 1.23, dtype=float)


def test_admet_predictor_predicts_with_mocked_models(monkeypatch):
    from amprenta_rag.ml.admet import predictor as predmod

    p = predmod.ADMETPredictor()

    # Avoid RDKit dependency
    monkeypatch.setattr(p, "_get_features", lambda smiles: np.zeros((2048 + 6,), dtype=np.float32))

    # Mock registry behavior
    dummy_cls_model = _DummyClassifier()
    dummy_reg_model = _DummyRegressor()

    dummy_record_cls = SimpleNamespace(id=uuid4())
    dummy_record_reg = SimpleNamespace(id=uuid4())

    def fake_get_active_model(name: str):  # noqa: ANN001
        if name == "admet_herg_xgb":
            return dummy_record_cls
        if name in ("admet_logs_xgb", "admet_logp_xgb"):
            return dummy_record_reg
        return None

    def fake_load_model(model_id):  # noqa: ANN001
        if model_id == dummy_record_cls.id:
            return dummy_cls_model
        return dummy_reg_model

    monkeypatch.setattr(p.registry, "get_active_model", fake_get_active_model)
    monkeypatch.setattr(p.registry, "load_model", fake_load_model)

    out = p.predict(["CCO"], endpoints=["herg", "logs"], include_shap=False)
    assert out and out[0]["smiles"] == "CCO"
    assert out[0]["herg"]["probability"] == pytest.approx(0.8)
    assert out[0]["logs"]["value"] == pytest.approx(1.23)


def test_admet_compound_endpoint_fetches_smiles_and_predicts(monkeypatch):
    from amprenta_rag.api.main import app
    import amprenta_rag.api.routers.chemistry as chem_router

    client = TestClient(app)
    compound_id = uuid4()

    # Mock DB session and compound lookup
    fake_db = SimpleNamespace(
        query=lambda model: SimpleNamespace(
            filter=lambda *args, **kwargs: SimpleNamespace(
                first=lambda: SimpleNamespace(id=compound_id, smiles="CCO")
            )
        )
    )

    @contextmanager
    def fake_db_session():
        yield fake_db

    monkeypatch.setattr(chem_router, "db_session", fake_db_session, raising=False)

    # Mock predictor
    class _Pred:
        def predict(self, smiles_list, endpoints=None, include_shap=False):  # noqa: ANN001
            return [{"smiles": smiles_list[0], "herg": {"probability": 0.5, "class": 0}}]

    monkeypatch.setattr(chem_router, "get_admet_predictor", lambda: _Pred(), raising=False)

    resp = client.post(
        f"/api/v1/chemistry/compounds/{compound_id}/predict-admet",
        json={"smiles_list": ["IGNORED"]},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["smiles"] == "CCO"
    assert body["herg"]["probability"] == 0.5


