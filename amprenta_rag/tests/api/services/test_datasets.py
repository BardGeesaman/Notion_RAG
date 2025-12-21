from __future__ import annotations

from enum import Enum
from types import SimpleNamespace
from uuid import uuid4

import pytest

from amprenta_rag.api.services import datasets


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)

    def any(self, **kwargs):
        return ("any", self.name, kwargs.get("id"))


class FakeDataset:
    id = Field("id")
    name = Field("name")
    omics_type = Field("omics_type")
    programs = Field("programs")
    experiments = Field("experiments")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.name = kwargs.get("name", "")
        self.omics_type = kwargs.get("omics_type", "")
        self.description = kwargs.get("description")
        self.file_paths = kwargs.get("file_paths", [])
        self.file_urls = kwargs.get("file_urls", [])
        self.organism = kwargs.get("organism", [])
        self.sample_type = kwargs.get("sample_type", [])
        self.disease = kwargs.get("disease", [])
        self.programs = kwargs.get("programs", [])
        self.experiments = kwargs.get("experiments", [])
        self.features = kwargs.get("features", [])


class FakeFeature:
    def __init__(self, fid, feature_type):
        self.id = fid
        self.feature_type = feature_type


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
        if model is FakeDataset:
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


class FakeUpdate:
    def __init__(self, **kwargs):
        self._data = kwargs

    def model_dump(self, exclude_unset: bool = False):
        return dict(self._data)


class FakeFeatureType(Enum):
    gene = "gene"
    protein = "protein"


@pytest.fixture(autouse=True)
def patch_model(monkeypatch):
    monkeypatch.setattr(datasets, "DatasetModel", FakeDataset)
    monkeypatch.setattr(datasets, "FeatureType", FakeFeatureType)
    yield


def test_create_dataset_adds_relationships(monkeypatch):
    db = FakeSession()
    program_id = uuid4()
    experiment_id = uuid4()
    fake_program = SimpleNamespace(id=program_id)
    fake_experiment = SimpleNamespace(id=experiment_id)
    monkeypatch.setattr(datasets, "uuid", SimpleNamespace(uuid4=lambda: uuid4()))
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: fake_program if pid == program_id else None)
    monkeypatch.setattr("amprenta_rag.api.services.experiments.get_experiment", lambda db, eid: fake_experiment if eid == experiment_id else None)

    payload = SimpleNamespace(
        name="D1",
        omics_type=SimpleNamespace(value="rna"),
        description="desc",
        file_paths=["p"],
        file_urls=["u"],
        organism=["human"],
        sample_type=["blood"],
        disease=["cancer"],
        program_ids=[program_id],
        experiment_ids=[experiment_id],
    )
    created = datasets.create_dataset(db, payload)

    assert created in db.items
    assert created.programs == [fake_program]
    assert created.experiments == [fake_experiment]
    assert created.name == "D1"
    assert created.omics_type == "rna"


def test_get_datasets_filters_and_limits():
    db = FakeSession()
    target_prog = SimpleNamespace(id=uuid4())
    other_prog = SimpleNamespace(id=uuid4())
    ds1 = FakeDataset(name="Alpha", omics_type="rna", programs=[target_prog])
    ds2 = FakeDataset(name="Beta", omics_type="dna", programs=[other_prog])
    db.items.extend([ds1, ds2])

    result = datasets.get_datasets(db, name_filter="alp", omics_type="rna", program_id=target_prog.id, limit=10)
    assert result == [ds1]

    with pytest.raises(ValueError):
        datasets.get_datasets(db, name_filter="x" * 101)


def test_update_dataset_updates_fields_and_relationships(monkeypatch):
    db = FakeSession()
    existing = FakeDataset(name="Old", omics_type="rna")
    db.items.append(existing)
    new_prog = SimpleNamespace(id=uuid4())
    new_exp = SimpleNamespace(id=uuid4())
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: new_prog)
    monkeypatch.setattr("amprenta_rag.api.services.experiments.get_experiment", lambda db, eid: new_exp)

    update = FakeUpdate(name="New", omics_type=SimpleNamespace(value="dna"), program_ids=[new_prog.id], experiment_ids=[new_exp.id])
    updated = datasets.update_dataset(db, existing.id, update)

    assert updated is existing
    assert existing.name == "New"
    assert existing.omics_type == "dna"
    assert existing.programs == [new_prog]
    assert existing.experiments == [new_exp]


def test_delete_dataset():
    db = FakeSession()
    existing = FakeDataset(name="Del")
    db.items.append(existing)

    assert datasets.delete_dataset(db, existing.id) is True
    assert existing not in db.items
    assert datasets.delete_dataset(db, uuid4()) is False


def test_get_dataset_features_by_type_groups_and_skips_invalid():
    db = FakeSession()
    f1 = FakeFeature(uuid4(), "gene")
    f2 = FakeFeature(uuid4(), "protein")
    f_invalid = FakeFeature(uuid4(), "unknown")
    ds = FakeDataset(features=[f1, f2, f_invalid])
    db.items.append(ds)

    grouped = datasets.get_dataset_features_by_type(db, ds.id)
    assert grouped[FakeFeatureType.gene] == [f1.id]
    assert grouped[FakeFeatureType.protein] == [f2.id]
    assert FakeFeatureType.gene in grouped and FakeFeatureType.protein in grouped

