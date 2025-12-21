from __future__ import annotations

from enum import Enum
from types import SimpleNamespace
from uuid import uuid4

import pytest

from amprenta_rag.api.services import signatures


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)

    def any(self, **kwargs):
        return ("any", self.name, kwargs.get("id"))


class FakeSignature:
    id = Field("id")
    name = Field("name")
    programs = Field("programs")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.name = kwargs.get("name", "")
        self.description = kwargs.get("description")
        self.components = []
        self.features = []
        self.programs = []
        self.modalities = []


class FakeSignatureComponent:
    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.signature_id = kwargs.get("signature_id")
        self.feature_id = kwargs.get("feature_id")
        self.feature_name = kwargs.get("feature_name")
        self.feature_type = kwargs.get("feature_type")
        self.direction = kwargs.get("direction")
        self.weight = kwargs.get("weight")


class FakeFeature:
    id = Field("id")

    def __init__(self, fid):
        self.id = fid


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
                    if kind == "eq" and actual != value:
                        ok = False
                        break
                    if kind == "ilike":
                        substring = value.replace("%", "").lower()
                        if substring not in (getattr(obj, field, "") or "").lower():
                            ok = False
                            break
                    if kind == "any":
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
        self.signatures = []
        self.features = []

    def query(self, model):
        if model is FakeSignature:
            return FakeQuery(self.signatures)
        if model is FakeFeature:
            return FakeQuery(self.features)
        return FakeQuery([])

    def add(self, obj):
        self.signatures.append(obj)

    def commit(self):
        return None

    def refresh(self, obj):
        return None

    def delete(self, obj):
        self.signatures = [s for s in self.signatures if s is not obj]


class FakeFeatureType(Enum):
    gene = "gene"
    protein = "protein"


class FakeDirection(Enum):
    up = "up"
    down = "down"


class FakeUpdate:
    def __init__(self, **kwargs):
        self._data = kwargs

    def model_dump(self, exclude_unset: bool = False):
        return dict(self._data)


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(signatures, "SignatureModel", FakeSignature)
    monkeypatch.setattr(signatures, "SignatureComponentModel", FakeSignatureComponent)
    monkeypatch.setattr(signatures, "FeatureModel", FakeFeature)
    yield


def _comp(feature_id=None, feature_name="f", ftype=FakeFeatureType.gene, feature_type=None, direction=None, weight=None, **_):
    ft = feature_type or ftype
    dir_obj = FakeDirection.up if direction == "up" else None
    return SimpleNamespace(
        feature_id=feature_id,
        feature_name=feature_name,
        feature_type=ft,
        direction=dir_obj,
        weight=weight,
    )


def test_create_signature_adds_components_features_and_programs(monkeypatch):
    db = FakeSession()
    feat_id = uuid4()
    db.features.append(FakeFeature(feat_id))
    program_id = uuid4()
    fake_program = SimpleNamespace(id=program_id)
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: fake_program if pid == program_id else None)
    payload = SimpleNamespace(
        name="Sig1",
        description="desc",
        components=[_comp(feature_id=feat_id, feature_type=FakeFeatureType.gene, weight=2.0)],
        program_ids=[program_id],
    )

    created = signatures.create_signature(db, payload)
    assert created in db.signatures
    assert created.components and created.components[0].feature_id == feat_id
    assert created.features == [db.features[0]]
    assert created.programs == [fake_program]
    assert "gene" in created.modalities


def test_get_signatures_filters():
    db = FakeSession()
    prog = SimpleNamespace(id=uuid4())
    s1 = FakeSignature(name="Alpha")
    s1.programs = [prog]
    s2 = FakeSignature(name="Beta")
    db.signatures.extend([s1, s2])

    result = signatures.get_signatures(db, name_filter="alp", program_id=prog.id)
    assert result == [s1]


def test_update_signature_replaces_components_and_programs(monkeypatch):
    db = FakeSession()
    existing = FakeSignature(name="Old")
    db.signatures.append(existing)
    feat_new = FakeFeature(uuid4())
    db.features.append(feat_new)
    new_prog = SimpleNamespace(id=uuid4())
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: new_prog)

    update = FakeUpdate(
        name="New",
        components=[_comp(feature_id=feat_new.id, feature_type=FakeFeatureType.protein, weight=3.0)],
        program_ids=[new_prog.id],
    )
    updated = signatures.update_signature(db, existing.id, update)

    assert updated is existing
    assert existing.name == "New"
    assert len(existing.components) == 1
    assert existing.features == [feat_new]
    assert existing.programs == [new_prog]
    assert set(existing.modalities) == {"protein"}


def test_delete_signature():
    db = FakeSession()
    sig = FakeSignature(name="Del")
    db.signatures.append(sig)

    assert signatures.delete_signature(db, sig.id) is True
    assert sig not in db.signatures
    assert signatures.delete_signature(db, uuid4()) is False

