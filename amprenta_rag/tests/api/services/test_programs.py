from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest

from amprenta_rag.api.services import programs


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)


class FakeProgram:
    id = Field("id")
    name = Field("name")
    description = Field("description")
    disease = Field("disease")

    def __init__(self, id, name, description, disease):
        self.id = id
        self.name = name
        self.description = description
        self.disease = disease


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
        if model is FakeProgram:
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


class FakeProgramUpdate:
    def __init__(self, **kwargs):
        self._data = kwargs

    def model_dump(self, exclude_unset: bool = False):
        return dict(self._data)


@pytest.fixture(autouse=True)
def patch_model(monkeypatch):
    monkeypatch.setattr(programs, "ProgramModel", FakeProgram)
    yield


def test_create_and_get_program():
    db = FakeSession()
    payload = SimpleNamespace(name="ProgA", description="desc", disease=["d1"])

    created = programs.create_program(db, payload)
    assert created.name == "ProgA"
    assert created.description == "desc"
    assert created.disease == ["d1"]
    assert created in db.items

    fetched = programs.get_program(db, created.id)
    assert fetched is created


def test_get_programs_with_filter():
    db = FakeSession()
    p1 = FakeProgram(uuid4(), "Alpha", "d", ["x"])
    p2 = FakeProgram(uuid4(), "Beta", "d", ["y"])
    db.items.extend([p1, p2])

    result = programs.get_programs(db, name_filter="alp")
    assert result == [p1]


def test_update_program_updates_fields_and_handles_none_list():
    db = FakeSession()
    existing = FakeProgram(uuid4(), "Old", "desc", ["a"])
    db.items.append(existing)

    update = FakeProgramUpdate(name="New", disease=None)
    updated = programs.update_program(db, existing.id, update)

    assert updated is existing
    assert existing.name == "New"
    assert existing.disease == []


def test_update_program_not_found_returns_none():
    db = FakeSession()
    update = FakeProgramUpdate(name="New")
    assert programs.update_program(db, uuid4(), update) is None


def test_delete_program():
    db = FakeSession()
    existing = FakeProgram(uuid4(), "Del", "desc", [])
    db.items.append(existing)

    assert programs.delete_program(db, existing.id) is True
    assert existing not in db.items
    assert programs.delete_program(db, uuid4()) is False

