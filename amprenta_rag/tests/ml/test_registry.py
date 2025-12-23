"""Unit tests for ML model registry."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict
from uuid import UUID, uuid4

import pytest


@dataclass
class DummyModel:
    """Pickle-friendly dummy model."""

    value: float = 1.0

    def predict(self, X):  # noqa: ANN001
        return [self.value for _ in range(len(X))]


class _FakeQuery:
    def __init__(self, store: Dict[UUID, Any], model_cls: Any):
        self._store = store
        self._model_cls = model_cls
        self._filters: Dict[str, Any] = {}

    def filter(self, *args, **kwargs):  # noqa: ANN001,ANN002,ANN003
        # Best-effort parse of SQLAlchemy BinaryExpressions: Model.col == value
        for expr in args:
            left = getattr(expr, "left", None)
            right = getattr(expr, "right", None)
            col = getattr(left, "name", None) if left is not None else None
            val = getattr(right, "value", None) if right is not None else None
            if col is not None:
                self._filters[col] = val
        self._filters.update(kwargs)
        return self

    def order_by(self, *args, **kwargs):  # noqa: ANN001,ANN002,ANN003
        return self

    def first(self):  # noqa: ANN001
        items = self.all()
        return items[0] if items else None

    def all(self):  # noqa: ANN001
        items = list(self._store.values())
        # Apply common filters used by MLModelRegistry
        if "id" in self._filters:
            items = [m for m in items if getattr(m, "id", None) == self._filters["id"]]
        if "name" in self._filters:
            items = [m for m in items if getattr(m, "name", None) == self._filters["name"]]
        if "status" in self._filters:
            items = [m for m in items if getattr(m, "status", None) == self._filters["status"]]
        if "model_type" in self._filters:
            items = [m for m in items if getattr(m, "model_type", None) == self._filters["model_type"]]
        # newest first if created_at exists
        items.sort(key=lambda m: getattr(m, "created_at", 0) or 0, reverse=True)
        return items

    def limit(self, n: int):  # noqa: ANN001
        self._limit = n
        return self


class _FakeDB:
    def __init__(self):
        self.store: Dict[UUID, Any] = {}

    def add(self, obj):  # noqa: ANN001
        # ID is usually assigned at DB time; set one for tests.
        if getattr(obj, "id", None) is None:
            obj.id = uuid4()
        self.store[obj.id] = obj

    def commit(self):  # noqa: D401
        """No-op commit."""

    def refresh(self, obj):  # noqa: ANN001
        # No-op refresh (object already in store)
        self.store[getattr(obj, "id")] = obj

    def query(self, model_cls):  # noqa: ANN001
        return _FakeQuery(self.store, model_cls)


@pytest.mark.unit
def test_registry_register_and_load(monkeypatch, tmp_path: Path):
    from amprenta_rag.ml.registry import MLModelRegistry
    from amprenta_rag import ml as _ml_pkg  # noqa: F401
    import amprenta_rag.ml.registry as regmod

    fake_db = _FakeDB()

    @contextmanager
    def fake_db_session():
        yield fake_db

    monkeypatch.setattr(regmod, "db_session", fake_db_session)

    reg = MLModelRegistry(artifact_base=tmp_path)
    model_obj = DummyModel(value=3.14)
    rec = reg.register_model(
        name="admet_herg_xgb",
        version="v1",
        model_type="admet_classification",
        framework="xgboost",
        model_object=model_obj,
        metrics={"auc": 0.9},
    )
    assert rec.id is not None

    loaded1 = reg.load_model(rec.id)
    loaded2 = reg.load_model(rec.id)
    assert isinstance(loaded1, DummyModel)
    assert loaded1.value == 3.14
    # Cached object identity
    assert loaded1 is loaded2


@pytest.mark.unit
def test_registry_list_models_filters(monkeypatch, tmp_path: Path):
    from amprenta_rag.ml.registry import MLModelRegistry
    import amprenta_rag.ml.registry as regmod

    fake_db = _FakeDB()

    @contextmanager
    def fake_db_session():
        yield fake_db

    monkeypatch.setattr(regmod, "db_session", fake_db_session)

    reg = MLModelRegistry(artifact_base=tmp_path)
    reg.register_model(
        name="m1",
        version="v1",
        model_type="admet_classification",
        framework="xgboost",
        model_object=DummyModel(1.0),
    )
    reg.register_model(
        name="m2",
        version="v1",
        model_type="admet_regression",
        framework="xgboost",
        model_object=DummyModel(2.0),
    )
    # archive one
    m3 = reg.register_model(
        name="m3",
        version="v1",
        model_type="admet_classification",
        framework="xgboost",
        model_object=DummyModel(3.0),
    )
    reg.archive_model(m3.id)

    all_models = reg.list_models()
    assert len(all_models) == 3

    active_cls = reg.list_models(model_type="admet_classification", status="active")
    assert {m.name for m in active_cls} == {"m1"}


def test_admet_predictor_features():
    """Test feature extraction from SMILES."""
    pytest.importorskip("rdkit")

    from amprenta_rag.ml.admet.predictor import ADMETPredictor

    predictor = ADMETPredictor()
    features = predictor._get_features("CCO")  # Ethanol
    assert features is not None
    assert len(features) == 2048 + 6  # Morgan FP + 6 descriptors


