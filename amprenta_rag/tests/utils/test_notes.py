from __future__ import annotations

from typing import Any, List
from uuid import UUID, uuid4

import pytest

from amprenta_rag.utils import notes


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def desc(self):
        return self


class FakeNote:
    entity_type = _Field("entity_type")
    entity_id = _Field("entity_id")
    id = _Field("id")
    created_at = _Field("created_at")

    def __init__(self, entity_type: str, entity_id: UUID, content: str, created_by_id: UUID | None):
        self._entity_type_val = entity_type
        self._entity_id_val = entity_id
        self.content = content
        self.created_by_id = created_by_id
        self.id: UUID | None = None
        self.created_at = None

    def __getattribute__(self, name: str):
        if name in ("entity_type", "entity_id"):
            return object.__getattribute__(self, f"_{name}_val")
        return object.__getattribute__(self, name)


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

    def order_by(self, *args: Any) -> "FakeQuery":
        return self

    def all(self) -> List[Any]:
        return list(self._data)

    def first(self) -> Any | None:
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self):
        self.notes: List[FakeNote] = []

    def add(self, obj: Any) -> None:
        if isinstance(obj, FakeNote):
            obj.id = obj.id or uuid4()
            self.notes.append(obj)

    def commit(self) -> None:
        return None

    def refresh(self, obj: Any) -> None:
        return None

    def query(self, model: Any) -> FakeQuery:
        if model is notes.Note:
            return FakeQuery(self.notes)
        return FakeQuery([])

    def delete(self, obj: Any) -> None:
        if isinstance(obj, FakeNote):
            self.notes = [n for n in self.notes if n is not obj]


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(notes, "Note", FakeNote)
    yield


def test_add_and_get_notes_happy_path():
    db = FakeSession()
    entity = uuid4()
    user = uuid4()

    note = notes.add_note("experiment", str(entity), "hello", str(user), db)
    assert note.entity_id == entity
    assert note.created_by_id == user

    retrieved = notes.get_notes("experiment", str(entity), db)
    assert len(retrieved) == 1
    assert retrieved[0].content == "hello"


def test_add_note_sets_created_by_optional():
    db = FakeSession()
    entity = uuid4()

    note = notes.add_note("compound", str(entity), "no user", "test", db)
    assert note.created_by_id is None


def test_delete_note_success_and_not_found():
    db = FakeSession()
    entity = uuid4()
    note = notes.add_note("signature", str(entity), "bye", str(uuid4()), db)

    assert notes.delete_note(str(note.id), db) is True
    assert notes.get_notes("signature", str(entity), db) == []

    assert notes.delete_note(str(uuid4()), db) is False

