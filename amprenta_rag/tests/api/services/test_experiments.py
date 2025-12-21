from __future__ import annotations

from types import SimpleNamespace
from uuid import uuid4

import pytest

from amprenta_rag.api.services import experiments


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def ilike(self, pattern: str):
        return ("ilike", self.name, pattern)

    def any(self, **kwargs):
        return ("any", self.name, kwargs.get("id"))


class FakeExperiment:
    id = Field("id")
    name = Field("name")
    programs = Field("programs")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.name = kwargs.get("name", "")
        self.type = kwargs.get("type")
        self.description = kwargs.get("description")
        self.disease = kwargs.get("disease", [])
        self.matrix = kwargs.get("matrix", [])
        self.model_systems = kwargs.get("model_systems", [])
        self.programs = kwargs.get("programs", [])


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
        if model is FakeExperiment:
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


@pytest.fixture(autouse=True)
def patch_model(monkeypatch):
    monkeypatch.setattr(experiments, "ExperimentModel", FakeExperiment)
    yield


def test_create_and_get_experiment(monkeypatch):
    db = FakeSession()
    program_id = uuid4()
    fake_program = SimpleNamespace(id=program_id)
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: fake_program if pid == program_id else None)
    payload = SimpleNamespace(
        name="E1",
        type="rna",
        description="desc",
        disease=["d"],
        matrix=["m"],
        model_systems=["s"],
        program_ids=[program_id],
    )

    created = experiments.create_experiment(db, payload)
    assert created.name == "E1"
    assert created.programs == [fake_program]
    assert experiments.get_experiment(db, created.id) is created


def test_get_experiments_filters():
    db = FakeSession()
    target_prog = SimpleNamespace(id=uuid4())
    exp1 = FakeExperiment(name="Alpha", programs=[target_prog])
    exp2 = FakeExperiment(name="Beta", programs=[])
    db.items.extend([exp1, exp2])

    result = experiments.get_experiments(db, name_filter="alp", program_id=target_prog.id)
    assert result == [exp1]


def test_update_experiment(monkeypatch):
    db = FakeSession()
    exp = FakeExperiment(name="Old", programs=[])
    db.items.append(exp)
    new_prog = SimpleNamespace(id=uuid4())
    monkeypatch.setattr("amprenta_rag.api.services.programs.get_program", lambda db, pid: new_prog)
    update = FakeUpdate(name="New", program_ids=[new_prog.id])

    updated = experiments.update_experiment(db, exp.id, update)
    assert updated is exp
    assert exp.name == "New"
    assert exp.programs == [new_prog]


def test_update_experiment_not_found():
    db = FakeSession()
    update = FakeUpdate(name="New")
    assert experiments.update_experiment(db, uuid4(), update) is None


def test_delete_experiment():
    db = FakeSession()
    exp = FakeExperiment(name="Del")
    db.items.append(exp)

    assert experiments.delete_experiment(db, exp.id) is True
    assert exp not in db.items
    assert experiments.delete_experiment(db, uuid4()) is False

