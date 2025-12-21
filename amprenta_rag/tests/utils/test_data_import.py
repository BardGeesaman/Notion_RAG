from __future__ import annotations

from typing import Any, List
from uuid import uuid4

import pandas as pd
import pytest

from amprenta_rag.utils import data_import


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def like(self, pattern: str):
        return ("like", self.name, pattern)

    def desc(self):
        return self


class FakeExperiment:
    id = _Field("id")
    name = _Field("name")

    def __init__(self, id: Any, name: str, **kwargs: Any):
        self._id_val = id
        self.name = name
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __getattribute__(self, name: str):
        if name == "id":
            return object.__getattribute__(self, "_id_val")
        return object.__getattribute__(self, name)


class FakeCompound:
    id = _Field("id")
    compound_id = _Field("compound_id")

    def __init__(self, compound_id: str, smiles: str = "CCO", inchi_key: str | None = None, **kwargs: Any):
        self.compound_id = compound_id
        self.smiles = smiles
        self.inchi_key = inchi_key
        self.canonical_smiles = smiles
        self.molecular_formula = "C2H6O"
        self.molecular_weight = None
        self.logp = None
        self.hbd_count = None
        self.hba_count = None
        self.rotatable_bonds = None
        self.external_ids = None
        self.created_at = None
        self.updated_at = None
        for k, v in kwargs.items():
            setattr(self, k, v)


class FakeSample:
    id = _Field("id")
    name = _Field("name")

    def __init__(self, name: str, **kwargs: Any):
        self.name = name
        self.id = uuid4()
        for k, v in kwargs.items():
            setattr(self, k, v)


class FakeStorageLocation:
    name = _Field("name")

    def __init__(self, name: str):
        self.name = name
        self.id = uuid4()


class FakeQuery:
    def __init__(self, data: List[Any]):
        self._data = list(data)
        self._order_desc = False

    def filter(self, *predicates: Any) -> "FakeQuery":
        data = list(self._data)
        for pred in predicates:
            if isinstance(pred, tuple):
                if pred[0] == "like":
                    _, field, pattern = pred
                    if pattern == "AMP-%":
                        data = [d for d in data if getattr(d, field).startswith("AMP-")]
                else:
                    field, val = pred
                    data = [d for d in data if getattr(d, field, None) == val]
        return FakeQuery(data)

    def order_by(self, *args: Any) -> "FakeQuery":
        return self

    def first(self) -> Any | None:
        return self._data[0] if self._data else None

    def all(self) -> List[Any]:
        return list(self._data)


class FakeSession:
    def __init__(self):
        self.experiments: List[FakeExperiment] = []
        self.compounds: List[FakeCompound] = []
        self.samples: List[FakeSample] = []
        self.storage_locations: List[FakeStorageLocation] = []
        self.added: List[Any] = []
        self.committed = False
        self.rolled_back = False

    def query(self, model: Any) -> FakeQuery:
        if model is data_import.Experiment:
            return FakeQuery(self.experiments)
        if model is data_import.Compound:
            return FakeQuery(self.compounds)
        if getattr(model, "__name__", "") == "StorageLocation":
            return FakeQuery(self.storage_locations)
        if model is data_import.Sample:
            return FakeQuery(self.samples)
        return FakeQuery([])

    def add(self, obj: Any) -> None:
        self.added.append(obj)
        if isinstance(obj, FakeExperiment):
            self.experiments.append(obj)
        if isinstance(obj, FakeCompound):
            self.compounds.append(obj)
        if isinstance(obj, FakeSample):
            self.samples.append(obj)

    def commit(self) -> None:
        self.committed = True

    def rollback(self) -> None:
        self.rolled_back = True


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(data_import, "Experiment", FakeExperiment)
    monkeypatch.setattr(data_import, "Compound", FakeCompound)
    monkeypatch.setattr(data_import, "Sample", FakeSample)
    # StorageLocation is imported inside function; patch module attribute for comparison
    monkeypatch.setattr(data_import, "StorageLocation", FakeStorageLocation, raising=False)
    yield


def test_validate_import_data_missing_columns():
    df = pd.DataFrame({"foo": [1]})
    errors = data_import.validate_import_data(df, "experiment")
    assert "Missing required columns" in errors[0]


def test_import_experiments_success():
    df = pd.DataFrame([{"name": "exp1", "design_metadata": None, "sample_groups": None}])
    db = FakeSession()
    result = data_import.import_experiments(df, db, user_id=str(uuid4()))
    assert result["created"] == 1
    assert result["errors"] == []
    assert isinstance(db.added[0], FakeExperiment)


def test_import_compounds_duplicate_skipped(monkeypatch):
    df = pd.DataFrame([{"smiles": "CCO"}])
    db = FakeSession()
    db.compounds.append(FakeCompound("AMP-00005"))

    monkeypatch.setattr(data_import, "check_duplicate", lambda s: True)
    result = data_import.import_compounds(df, db)
    assert result["created"] == 0
    assert result["duplicates"] == 1


def test_import_compounds_generates_id_and_descriptors(monkeypatch):
    df = pd.DataFrame([{"smiles": "CCO"}])
    db = FakeSession()
    db.compounds.append(FakeCompound("AMP-00001"))

    monkeypatch.setattr(data_import, "check_duplicate", lambda s: False)
    monkeypatch.setattr(data_import, "normalize_smiles", lambda s: ("CCO", "INCHI", "C2H6O"))
    monkeypatch.setattr(
        data_import,
        "compute_molecular_descriptors",
        lambda s: {"molecular_weight": 46.0, "logp": 0.1, "hbd_count": 1, "hba_count": 1, "rotatable_bonds": 0},
    )

    result = data_import.import_compounds(df, db)
    assert result["created"] == 1
    assert result["duplicates"] == 0
    added = db.added[-1]
    assert added.compound_id.startswith("AMP-")
    assert added.molecular_weight == 46.0


def test_import_samples_missing_refs_records_errors(monkeypatch):
    df = pd.DataFrame([{"name": "sample1", "storage_location": "L1", "experiment": "E1", "parent_sample": "P1"}])
    db = FakeSession()
    # reference objects not present
    result = data_import.import_samples(df, db, user_id=str(uuid4()))
    assert result["created"] == 1
    assert len(result["errors"]) == 3  # missing storage, experiment, parent


def test_get_template_returns_expected_columns():
    tmpl = data_import.get_template("compound")
    assert list(tmpl.columns)[:3] == ["compound_id", "smiles", "inchi_key"]

