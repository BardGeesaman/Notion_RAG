from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List
from uuid import UUID, uuid4

import pytest

from amprenta_rag.utils import comments


class _Field:
    """Lightweight stand-in for SQLAlchemy column attributes."""

    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):  # type: ignore[override]
        return (self.name, other)

    def is_(self, other: object):
        return ("is", self.name, other)

    def desc(self) -> "_Field":
        return self


class FakeComment:
    entity_type = _Field("entity_type")
    entity_id = _Field("entity_id")
    parent_id = _Field("parent_id")
    created_at = _Field("created_at")
    created_by_id = _Field("created_by_id")
    id = _Field("id")

    def __init__(
        self,
        entity_type: str,
        entity_id: UUID,
        content: str,
        created_by_id: UUID | None,
        parent_id: UUID | None = None,
        created_at: Any = None,
    ):
        self._entity_type_val = entity_type
        self._entity_id_val = entity_id
        self._content_val = content
        self._created_by_id_val = created_by_id
        self._parent_id_val = parent_id
        self._created_at_val = created_at
        self.id: UUID | None = None

    def __getattribute__(self, name: str):
        if name in ("entity_type", "entity_id", "parent_id", "created_at", "created_by_id"):
            return object.__getattribute__(self, f"_{name}_val")
        if name == "content":
            return object.__getattribute__(self, "_content_val")
        return object.__getattribute__(self, name)


class FakeUser:
    id = _Field("id")

    def __init__(self, id: UUID, username: str):
        self._id_val = id
        self.username = username

    def __getattribute__(self, name: str):
        if name == "id":
            return object.__getattribute__(self, "_id_val")
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
                    if len(pred) == 3 and pred[0] == "is":
                        _, field, val = pred
                        if getattr(obj, field, None) is not val:
                            match = False
                            break
                    else:
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

    def all(self) -> List[Any]:
        return list(self._data)

    def first(self) -> Any | None:
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self):
        self.comments: List[FakeComment] = []
        self.users: List[FakeUser] = []

    def add(self, obj: Any) -> None:
        if isinstance(obj, FakeComment):
            obj.id = obj.id or uuid4()
            self.comments.append(obj)
        elif isinstance(obj, FakeUser):
            self.users.append(obj)

    def commit(self) -> None:
        return None

    def refresh(self, obj: Any) -> None:
        return None

    def query(self, model: Any) -> FakeQuery:
        if model is comments.Comment:
            return FakeQuery(self.comments)
        if model is comments.User:
            return FakeQuery(self.users)
        return FakeQuery([])

    def delete(self, obj: Any) -> None:
        if isinstance(obj, FakeComment):
            self.comments = [c for c in self.comments if c is not obj]


@pytest.fixture(autouse=True)
def patch_models(monkeypatch):
    monkeypatch.setattr(comments, "Comment", FakeComment)
    monkeypatch.setattr(comments, "User", FakeUser)
    yield


def test_add_and_get_comments_with_replies():
    db = FakeSession()
    user = FakeUser(id=uuid4(), username="alice")
    db.add(user)
    entity = uuid4()

    parent = comments.add_comment("experiment", entity, "root", user.id, db)
    reply = comments.add_comment("experiment", entity, "child", user.id, db, parent_id=parent.id)

    results = comments.get_comments("experiment", entity, db)
    assert len(results) == 1
    top = results[0]
    assert top["id"] == str(parent.id)
    assert top["author"] == "alice"
    assert len(top["replies"]) == 1
    assert top["replies"][0]["id"] == str(reply.id)
    assert top["replies"][0]["author"] == "alice"


def test_get_comments_unknown_user_and_order():
    db = FakeSession()
    entity = uuid4()
    # No user added for this comment
    comments.add_comment("signature", entity, "orphan", uuid4(), db)

    results = comments.get_comments("signature", entity, db)
    assert results[0]["author"] == "Unknown"
    assert results[0]["replies"] == []


def test_delete_comment_only_by_author():
    db = FakeSession()
    user_id = uuid4()
    other_id = uuid4()
    entity = uuid4()

    parent = comments.add_comment("dataset", entity, "to-delete", user_id, db)
    assert comments.delete_comment(parent.id, user_id, db) is True
    assert comments.get_comments("dataset", entity, db) == []

    # Recreate and try with different user
    parent = comments.add_comment("dataset", entity, "keep", user_id, db)
    assert comments.delete_comment(parent.id, other_id, db) is False
    assert comments.get_comments("dataset", entity, db)[0]["id"] == str(parent.id)

