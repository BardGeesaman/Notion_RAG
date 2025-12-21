from __future__ import annotations

from enum import Enum
from types import SimpleNamespace
from uuid import uuid4

import pytest

from amprenta_rag.api.services import features


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)

    def any(self, **kwargs):
        return ("any", self.name, kwargs.get("id"))


class FakeFeature:
    id = Field("id")
    name = Field("name")
    feature_type = Field("feature_type")
    datasets = Field("datasets")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.name = kwargs.get("name", "")
        self.feature_type = kwargs.get("feature_type", "")
        self.normalized_name = kwargs.get("normalized_name")
        self.aliases = kwargs.get("aliases", [])
        self.external_ids = kwargs.get("external_ids", {})
        self.datasets = kwargs.get("datasets", [])


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, *predicates):
        filtered = []
        for obj in self._data:
            ok = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    kind, field, value = pred
                    actual = getattr(obj, field, None)
                    if kind == "eq":
                        if actual != value:
                            ok = False
                            break
                    elif kind == "ilike":
                        substring = value.replace("%", "").lower()
                        if substring not in (getattr(obj, field, "") or "").lower():
                            ok = False
                            break
                    elif kind == "any":
                        coll = getattr(obj, field, [])
                        if not any(getattr(item, "id", None) == value for item in coll):
                            ok = False
                            break
                else:
                    ok = False
                    break
            if ok:
                filtered.append(obj)
        return FakeQuery(filtered)

    def offset(self, n: int):
        return FakeQuery(self._data[n:])

    def limit(self, n: int):
        return FakeQuery(self._data[:n])

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self):
        self.items = []
        self.deleted = []

    def query(self, model):
        if model is FakeFeature:
            return FakeQuery(self.items)
        return FakeQuery([])

    def add(self, obj):
        self.items.append(obj)

    def commit(self):
        return None

    def refresh(self, obj):
        return None

    def delete(self, obj):
        self.deleted.append(obj)
        self.items = [i for i in self.items if i is not obj]


class FakeFeatureType(Enum):
    gene = "gene"
    protein = "protein"


class FakeUpdate:
    def __init__(self, **kwargs):
        self._data = kwargs

    def model_dump(self, exclude_unset: bool = False):
        return dict(self._data)


@pytest.fixture(autouse=True)
def patch_model(monkeypatch):
    monkeypatch.setattr(features, "FeatureModel", FakeFeature)
    monkeypatch.setattr(features, "FeatureType", FakeFeatureType)
    yield


def test_create_and_get_feature():
    db = FakeSession()
    payload = SimpleNamespace(
        name="F1",
        feature_type=SimpleNamespace(value="gene"),
        normalized_name="nF1",
        aliases=["a"],
        external_ids={"k": "v"},
    )
    created = features.create_feature(db, payload)
    assert created.name == "F1"
    assert created.feature_type == "gene"
    assert features.get_feature(db, created.id) is created


def test_get_feature_by_name_with_type_filter():
    db = FakeSession()
    f1 = FakeFeature(name="F1", feature_type="gene")
    f2 = FakeFeature(name="F1", feature_type="protein")
    db.items.extend([f1, f2])

    result = features.get_feature_by_name(db, "F1", feature_type=FakeFeatureType.gene)
    assert result is f1
    assert features.get_feature_by_name(db, "missing") is None


def test_get_features_filters():
    db = FakeSession()
    ds = SimpleNamespace(id=uuid4())
    f1 = FakeFeature(name="Alpha", feature_type="gene", datasets=[ds])
    f2 = FakeFeature(name="Beta", feature_type="protein", datasets=[])
    db.items.extend([f1, f2])

    filtered = features.get_features(db, name_filter="alp", feature_type=FakeFeatureType.gene, dataset_id=ds.id)
    assert filtered == [f1]


def test_update_feature_and_delete():
    db = FakeSession()
    f1 = FakeFeature(name="Old", feature_type="gene")
    db.items.append(f1)
    update = FakeUpdate(name="New", aliases=["a"])

    updated = features.update_feature(db, f1.id, update)
    assert updated is f1
    assert f1.name == "New"
    assert f1.aliases == ["a"]

    assert features.delete_feature(db, f1.id) is True
    assert f1 not in db.items
    assert features.delete_feature(db, uuid4()) is False

