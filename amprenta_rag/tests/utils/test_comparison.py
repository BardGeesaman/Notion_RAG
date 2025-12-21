from __future__ import annotations

from uuid import uuid4

import pytest

from amprenta_rag.utils import comparison


class FakeExperiment:
    def __init__(self, id, name, description, design_type):
        self.id = id
        self.name = name
        self.description = description
        self.design_type = design_type


class FakeCompound:
    def __init__(self, id, compound_id, smiles, mw, logp):
        self.id = id
        self.compound_id = compound_id
        self.smiles = smiles
        self.molecular_weight = mw
        self.logp = logp


class FakeQuery:
    def __init__(self, obj):
        self._obj = obj

    def filter(self, *args, **kwargs):
        return self

    def first(self):
        return self._obj


class FakeSession:
    def __init__(self, obj1, obj2):
        self.obj1 = obj1
        self.obj2 = obj2
        self.call = 0

    def query(self, model):
        self.call += 1
        return FakeQuery(self.obj1 if self.call == 1 else self.obj2)


def test_compare_experiments_differences():
    exp1 = FakeExperiment(uuid4(), "a", "d1", "case_control")
    exp2 = FakeExperiment(uuid4(), "b", "d2", "time_course")
    db = FakeSession(exp1, exp2)
    result = comparison.compare_experiments(str(exp1.id), str(exp2.id), db)
    assert set(result["differences"]) == {"name", "description", "design_type"}
    assert result["item1"]["name"] == "a"


def test_compare_compounds_matches():
    comp1 = FakeCompound(uuid4(), "C1", "CCO", 46.0, 0.1)
    comp2 = FakeCompound(uuid4(), "C1", "CCO", 46.0, 0.1)
    db = FakeSession(comp1, comp2)
    result = comparison.compare_compounds(str(comp1.id), str(comp2.id), db)
    assert result["differences"] == []


def test_compare_compounds_raises_on_missing():
    comp1 = None
    comp2 = FakeCompound(uuid4(), "C1", "CCO", 46.0, 0.1)
    db = FakeSession(comp1, comp2)
    with pytest.raises(ValueError):
        comparison.compare_compounds(str(uuid4()), str(uuid4()), db)

