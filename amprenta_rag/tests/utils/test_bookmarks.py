from __future__ import annotations

from typing import List, Any
from uuid import UUID, uuid4

import pytest

from amprenta_rag.utils import bookmarks


class _Field:
    """Lightweight stand-in for SQLAlchemy column attributes."""

    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def desc(self) -> "_Field":
        return self


class FakeBookmark:
    """Simple bookmark record used for tests."""

    # class-level fields to support comparisons in filters/order_by
    user_id = _Field("user_id")
    entity_type = _Field("entity_type")
    entity_id = _Field("entity_id")
    created_at = _Field("created_at")

    def __init__(self, user_id: UUID, entity_type: str, entity_id: UUID, id: UUID | None = None):
        self._user_id_val = user_id
        self._entity_type_val = entity_type
        self._entity_id_val = entity_id
        self.id: UUID | None = id

    def __getattribute__(self, name: str):
        if name in ("user_id", "entity_type", "entity_id"):
            return object.__getattribute__(self, f"_{name}_val")
        return object.__getattribute__(self, name)


class FakeQuery:
    def __init__(self, data: List[FakeBookmark]):
        self._data = data

    def filter(self, *predicates: Any) -> "FakeQuery":
        filtered: List[FakeBookmark] = []
        for obj in self._data:
            match = True
            for pred in predicates:
                if isinstance(pred, tuple):
                    field, val = pred
                    if getattr(obj, field, None) != val:
                        match = False
                        break
                elif pred is False:
                    match = False
                    break
            if match:
                filtered.append(obj)
        return FakeQuery(filtered)

    def order_by(self, *args: Any) -> "FakeQuery":
        return self

    def first(self) -> FakeBookmark | None:
        return self._data[0] if self._data else None

    def all(self) -> List[FakeBookmark]:
        return list(self._data)


class FakeSession:
    def __init__(self):
        self.data: List[FakeBookmark] = []

    def add(self, obj: FakeBookmark) -> None:
        obj.id = obj.id or uuid4()
        self.data.append(obj)

    def commit(self) -> None:
        return None

    def refresh(self, obj: FakeBookmark) -> None:
        return None

    def query(self, model: Any) -> FakeQuery:
        # All bookmark operations target Bookmark model; return a query over stored items.
        return FakeQuery(self.data)

    def delete(self, obj: FakeBookmark) -> None:
        self.data = [b for b in self.data if b is not obj]


@pytest.fixture(autouse=True)
def patch_bookmark(monkeypatch):
    monkeypatch.setattr(bookmarks, "Bookmark", FakeBookmark)
    yield


def test_add_and_get_bookmark_happy_path():
    db = FakeSession()
    user_id = uuid4()
    entity_id = uuid4()

    bm = bookmarks.add_bookmark(str(user_id), "experiment", str(entity_id), db)
    assert bm.user_id == user_id
    assert bm.entity_type == "experiment"
    assert bm.entity_id == entity_id
    assert bm.id is not None

    result = bookmarks.get_user_bookmarks(str(user_id), db)
    assert len(result) == 1
    stored = result[0]
    assert stored.user_id == user_id
    assert stored.entity_type == "experiment"
    assert stored.entity_id == entity_id


def test_remove_bookmark():
    db = FakeSession()
    user_id = uuid4()
    entity_id = uuid4()
    bookmarks.add_bookmark(user_id, "signature", entity_id, db)
    assert bookmarks.is_bookmarked(user_id, "signature", entity_id, db)

    bookmarks.remove_bookmark(user_id, "signature", entity_id, db)
    assert not bookmarks.is_bookmarked(user_id, "signature", entity_id, db)
    assert bookmarks.get_user_bookmarks(user_id, db) == []


def test_is_bookmarked_false_when_missing():
    db = FakeSession()
    user_id = uuid4()
    entity_id = uuid4()
    assert not bookmarks.is_bookmarked(user_id, "dataset", entity_id, db)


def test_add_bookmark_accepts_uuid_inputs():
    db = FakeSession()
    user_id = uuid4()
    entity_id = uuid4()
    bookmarks.add_bookmark(user_id, "dataset", entity_id, db)
    assert bookmarks.is_bookmarked(user_id, "dataset", entity_id, db)


def test_get_user_bookmarks_empty():
    db = FakeSession()
    user_id = uuid4()
    assert bookmarks.get_user_bookmarks(user_id, db) == []

