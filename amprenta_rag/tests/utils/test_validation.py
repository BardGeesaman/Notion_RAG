from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List
from uuid import uuid4

import pytest

from amprenta_rag.utils import validation


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)


@dataclass
class FakeExperiment:
    id: Any
    name: str | None = None
    design_type: str | None = None


@dataclass
class FakeCompound:
    id: Any
    smiles: str | None = None
    molecular_weight: float | None = None


class FakeQuery:
    def __init__(self, data: List[Any]):
        self._data = list(data)

    def filter(self, *predicates: Any) -> "FakeQuery":
        filtered: List[Any] = []
        for obj in self._data:
            match = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    field, val = pred
                    if getattr(obj, field, None) != val:
                        match = False
                        break
            if match:
                filtered.append(obj)
        return FakeQuery(filtered)

    def all(self) -> List[Any]:
        return list(self._data)

    def first(self) -> Any | None:
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self):
        self.experiments: List[FakeExperiment] = []
        self.compounds: List[FakeCompound] = []

    def query(self, model: Any) -> FakeQuery:
        if model is validation.Experiment:
            return FakeQuery(self.experiments)
        if model is validation.Compound:
            return FakeQuery(self.compounds)
        return FakeQuery([])


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(validation, "Experiment", FakeExperiment)
    monkeypatch.setattr(validation, "Compound", FakeCompound)
    yield


def test_validate_experiment_not_found():
    db = FakeSession()
    issues = validation.validate_experiment(str(uuid4()), db)
    assert len(issues) == 1
    assert issues[0].severity == "error"
    assert issues[0].issue == "Experiment not found"


def test_validate_experiment_name_and_design_type():
    db = FakeSession()
    bad = FakeExperiment(id=uuid4(), name="  ", design_type="invalid")
    db.experiments.append(bad)
    issues = validation.validate_experiment(str(bad.id), db)
    assert {i.field for i in issues} == {"name", "design_type"}
    assert any("Invalid design_type" in i.issue for i in issues)


def test_validate_compound_not_found():
    db = FakeSession()
    issues = validation.validate_compound(str(uuid4()), db)
    assert len(issues) == 1
    assert issues[0].issue == "Compound not found"


def test_validate_compound_missing_smiles_and_weight_warning(monkeypatch):
    db = FakeSession()
    comp = FakeCompound(id=uuid4(), smiles=None, molecular_weight=None)
    db.compounds.append(comp)
    issues = validation.validate_compound(str(comp.id), db)
    assert any(i.field == "smiles" and i.severity == "error" for i in issues)
    assert any(i.field == "molecular_weight" for i in issues) is False


def test_validate_compound_invalid_smiles(monkeypatch):
    db = FakeSession()
    comp = FakeCompound(id=uuid4(), smiles="not-a-smiles", molecular_weight=None)
    db.compounds.append(comp)

    class FakeChem:
        @staticmethod
        def MolFromSmiles(smiles: str):
            return None

    monkeypatch.setitem(validation.__dict__, "Chem", FakeChem)
    issues = validation.validate_compound(str(comp.id), db)
    assert any(i.field == "smiles" and "Invalid SMILES" in i.issue for i in issues)


def test_validate_compound_molecular_weight_warning(monkeypatch):
    db = FakeSession()
    comp = FakeCompound(id=uuid4(), smiles="CCO", molecular_weight=None)
    db.compounds.append(comp)

    class FakeChem:
        @staticmethod
        def MolFromSmiles(smiles: str):
            return object()

    monkeypatch.setitem(validation.__dict__, "Chem", FakeChem)
    issues = validation.validate_compound(str(comp.id), db)
    assert any(i.field == "molecular_weight" and i.severity == "warning" for i in issues)


def test_run_all_validations_aggregates(monkeypatch):
    db = FakeSession()
    exp = FakeExperiment(id=uuid4(), name=None, design_type=None)
    comp = FakeCompound(id=uuid4(), smiles=None, molecular_weight=None)
    db.experiments.append(exp)
    db.compounds.append(comp)

    # Avoid RDKit import attempt
    monkeypatch.setitem(validation.__dict__, "Chem", type("FakeChem", (), {"MolFromSmiles": staticmethod(lambda _: None)}))

    issues = validation.run_all_validations(db)
    assert len(issues) >= 2
    fields = {i.field for i in issues}
    assert "name" in fields or "smiles" in fields

