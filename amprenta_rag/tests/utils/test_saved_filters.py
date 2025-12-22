from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List
from uuid import UUID, uuid4

import pytest

from amprenta_rag.utils import saved_filters


class _Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def desc(self):
        return self


class FakeSavedFilter:
    # class-level field stand-ins for query building
    entity_type = _Field("entity_type")
    user_id = _Field("user_id")
    id = _Field("id")
    created_at = _Field("created_at")

    def __init__(self, name: str, entity_type: str, filters: Dict[str, Any], user_id: UUID, id: UUID | None = None):
        self.name = name
        self._entity_type_val = entity_type
        self._filters_val = filters
        self._user_id_val = user_id
        self.filters = filters
        self._id_val = id

    def __getattribute__(self, name: str):
        if name in ("entity_type", "user_id", "id"):
            stored = object.__getattribute__(self, "__dict__").get(f"_{name}_val", None)
            if stored is not None:
                return stored
            return object.__getattribute__(self, name)
        return object.__getattribute__(self, name)


class FakeQuery:
    def __init__(self, data: List[FakeSavedFilter]):
        self._data = data

    def filter(self, *preds: Any) -> "FakeQuery":
        filtered: List[FakeSavedFilter] = []
        for obj in self._data:
            match = True
            for p in preds:
                if isinstance(p, tuple):
                    field, val = p
                    # Use stored values for comparisons
                    actual = getattr(obj, f"_{field}_val", getattr(obj, field, None))
                    if actual != val:
                        match = False
                        break
            if match:
                filtered.append(obj)
        return FakeQuery(filtered)

    def order_by(self, *args: Any) -> "FakeQuery":
        return self

    def first(self) -> FakeSavedFilter | None:
        return self._data[0] if self._data else None

    def all(self) -> List[FakeSavedFilter]:
        return list(self._data)


class FakeSession:
    def __init__(self):
        self.data: List[FakeSavedFilter] = []

    def add(self, obj: FakeSavedFilter) -> None:
        obj._id_val = obj._id_val or uuid4()
        self.data.append(obj)

    def commit(self) -> None:
        return None

    def refresh(self, obj: FakeSavedFilter) -> None:
        return None

    def query(self, model: Any) -> FakeQuery:
        return FakeQuery(self.data)

    def delete(self, obj: FakeSavedFilter) -> None:
        self.data = [f for f in self.data if f is not obj]


@pytest.fixture(autouse=True)
def patch_model(monkeypatch):
    monkeypatch.setattr(saved_filters, "SavedFilter", FakeSavedFilter)
    yield


def test_save_and_get_filters():
    db = FakeSession()
    user_id = uuid4()
    sf = saved_filters.save_filter("Name", "experiment", {"k": "v"}, str(user_id), db)
    assert sf.name == "Name"
    assert sf.entity_type == "experiment"
    assert sf.filters == {"k": "v"}
    assert sf.user_id == user_id

    fetched = saved_filters.get_user_filters("experiment", str(user_id), db)
    assert len(fetched) == 1
    assert fetched[0].name == "Name"


def test_get_user_filters_empty():
    db = FakeSession()
    user_id = uuid4()
    fetched = saved_filters.get_user_filters("compound", str(user_id), db)
    assert fetched == []


def test_delete_filter_success():
    db = FakeSession()
    user_id = uuid4()
    sf = saved_filters.save_filter("F1", "dataset", {"a": 1}, user_id, db)
    assert saved_filters.delete_filter(str(sf.id), db) is True
    assert saved_filters.get_user_filters("dataset", user_id, db) == []


def test_delete_filter_not_found():
    db = FakeSession()
    assert saved_filters.delete_filter(str(uuid4()), db) is False

