from __future__ import annotations

from uuid import uuid4

from amprenta_rag.utils import global_search


class _Col:
    def __init__(self, name: str):
        self.name = name

    def ilike(self, val):
        return ("ilike", self.name, val)

    def __or__(self, other):
        return ("or", self, other)


class FakeExperiment:
    id = uuid4()
    name = _Col("name")
    description = _Col("description")

    def __init__(self, name: str, description: str = "", type_: str = "t"):
        self.id = uuid4()
        self.name = name
        self.description = description
        self.type = type_


class FakeCompound:
    id = uuid4()
    compound_id = _Col("compound_id")
    canonical_smiles = _Col("canonical_smiles")

    def __init__(self, cid: str, smiles: str):
        self.id = uuid4()
        self.compound_id = cid
        self.canonical_smiles = smiles
        self.smiles = smiles


class FakeSignature:
    id = uuid4()
    name = _Col("name")

    def __init__(self, name: str):
        self.id = uuid4()
        self.name = name


class FakeDataset:
    id = uuid4()
    name = _Col("name")

    def __init__(self, name: str):
        self.id = uuid4()
        self.name = name


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, *args, **kwargs):
        return self

    def limit(self, n: int):
        return self

    def all(self):
        return list(self._data)


class FakeSession:
    def __init__(self, experiments=None, compounds=None, signatures=None, datasets=None):
        self.experiments = experiments or []
        self.compounds = compounds or []
        self.signatures = signatures or []
        self.datasets = datasets or []

    def query(self, model):
        if model is global_search.Experiment:
            return FakeQuery(self.experiments)
        if model is global_search.Compound:
            return FakeQuery(self.compounds)
        if model is global_search.Signature:
            return FakeQuery(self.signatures)
        if model is global_search.Dataset:
            return FakeQuery(self.datasets)
        return FakeQuery([])


def test_global_search_returns_grouped_results(monkeypatch):
    monkeypatch.setattr(global_search, "Experiment", FakeExperiment)
    monkeypatch.setattr(global_search, "Compound", FakeCompound)
    monkeypatch.setattr(global_search, "Signature", FakeSignature)
    monkeypatch.setattr(global_search, "Dataset", FakeDataset)

    db = FakeSession(
        experiments=[FakeExperiment("exp apple", "desc")],
        compounds=[FakeCompound("C1", "apple")],
        signatures=[FakeSignature("sig apple")],
        datasets=[FakeDataset("ds apple")],
    )

    results = global_search.global_search("apple", db, limit=2)
    assert len(results["experiments"]) == 1
    assert len(results["compounds"]) == 1
    assert len(results["signatures"]) == 1
    assert len(results["datasets"]) == 1


def test_global_search_empty_query_returns_empty(monkeypatch):
    db = FakeSession()
    out = global_search.global_search("   ", db)
    assert all(not v for v in out.values())

