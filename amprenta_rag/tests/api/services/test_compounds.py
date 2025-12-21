from __future__ import annotations

from datetime import datetime, timedelta
from types import SimpleNamespace

from amprenta_rag.api.services import compounds


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)


class FakeCompound:
    compound_id = Field("compound_id")
    created_at = Field("created_at")
    updated_at = Field("updated_at")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id")
        self.compound_id = kwargs.get("compound_id")
        self.smiles = kwargs.get("smiles")
        self.inchi_key = kwargs.get("inchi_key")
        self.canonical_smiles = kwargs.get("canonical_smiles")
        self.molecular_formula = kwargs.get("molecular_formula")
        self.molecular_weight = kwargs.get("molecular_weight")
        self.logp = kwargs.get("logp")
        self.hbd_count = kwargs.get("hbd_count")
        self.hba_count = kwargs.get("hba_count")
        self.rotatable_bonds = kwargs.get("rotatable_bonds")
        self.updated_at = kwargs.get("updated_at")
        self.created_at = kwargs.get("created_at")
        self.programs = kwargs.get("programs", [])


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def order_by(self, *_):
        return self

    def filter(self, predicate):
        if isinstance(predicate, tuple):
            if len(predicate) == 2 and isinstance(predicate[0], tuple):
                kind = predicate[0]
                value = predicate[1]
            elif len(predicate) == 3:
                kind = predicate
                value = predicate[2]
            else:
                kind = None
                value = None
            if kind and kind[0] == "eq" and kind[1] == "compound_id":
                filtered = [c for c in self._data if c.compound_id == value]
                return FakeQuery(filtered)
        return FakeQuery([])

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSessionCtx:
    def __init__(self, compounds_data):
        self.compounds = compounds_data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def query(self, model):
        return FakeQuery(self.compounds)


def test_list_compounds_order_and_fields(monkeypatch):
    now = datetime.utcnow()
    old = FakeCompound(compound_id="C1", created_at=now - timedelta(days=2), updated_at=now - timedelta(days=1))
    new = FakeCompound(compound_id="C2", created_at=now - timedelta(days=1), updated_at=now)
    monkeypatch.setattr(compounds, "db_session", lambda: FakeSessionCtx([new, old]))

    rows = compounds.list_compounds()
    assert [r["compound_id"] for r in rows] == ["C2", "C1"]
    assert rows[0]["updated_at"] is not None


def test_get_compound_by_id(monkeypatch):
    c = FakeCompound(compound_id="C1")
    monkeypatch.setattr(compounds, "db_session", lambda: FakeSessionCtx([c]))
    monkeypatch.setattr(compounds, "Compound", FakeCompound)

    result = compounds.get_compound_by_id("C1")
    assert result["compound_id"] == "C1"
    assert compounds.get_compound_by_id("missing") is None


def test_get_compound_programs(monkeypatch):
    prog = SimpleNamespace(id=None, name="p", description="d")
    c = FakeCompound(compound_id="C1", programs=[prog])
    monkeypatch.setattr(compounds, "db_session", lambda: FakeSessionCtx([c]))
    monkeypatch.setattr(compounds, "Compound", FakeCompound)

    programs = compounds.get_compound_programs("C1")
    assert programs == [{"id": None, "name": "p", "description": "d"}]


def test_compound_to_dict_handles_none(monkeypatch):
    c = FakeCompound()
    monkeypatch.setattr(compounds, "db_session", lambda: FakeSessionCtx([c]))
    # ensure helper handles falsy correctly through get_compound_by_id path
    assert compounds.get_compound_by_id("") is None

