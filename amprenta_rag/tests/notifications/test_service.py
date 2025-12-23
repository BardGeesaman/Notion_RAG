from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any, List, Optional
from uuid import UUID, uuid4


from amprenta_rag.notifications import service as ns


class _Col:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: object):
        return (self.name, "==", other)

    def is_(self, other: object):
        return (self.name, "is", other)

    def desc(self):
        return (self.name, "desc")


class _Query:
    def __init__(self, db: "_DB"):
        self._db = db
        self._filters: list[Any] = []
        self._limit: Optional[int] = None
        self._order_desc: bool = False

    def filter(self, *args: Any):
        self._filters.extend(args)
        self._db.last_filters = list(args)
        return self

    def order_by(self, *_args: Any):
        self._order_desc = True
        return self

    def limit(self, n: int):
        self._limit = n
        return self

    def all(self) -> List[Any]:
        # Only used by get_unread_notifications
        rows = [n for n in self._db.notifications if (not n.is_read)]
        if self._order_desc:
            rows.sort(key=lambda x: x.created_at, reverse=True)
        if self._limit is not None:
            rows = rows[: self._limit]
        return rows

    def count(self) -> int:
        return len([n for n in self._db.notifications if not n.is_read])

    def first(self):
        # Only used by mark_as_read: match by id if present
        for arg in self._filters:
            if isinstance(arg, tuple) and arg[0] == "id" and arg[1] == "==" and isinstance(arg[2], UUID):
                return self._db.by_id.get(arg[2])
        return None

    def update(self, values: dict[str, Any]) -> int:
        # Only used by mark_all_read
        updated = 0
        for n in self._db.notifications:
            if not n.is_read:
                n.is_read = bool(values.get("is_read", n.is_read))
                updated += 1
        return updated


class _DB:
    def __init__(self):
        self.notifications: list[Any] = []
        self.by_id: dict[UUID, Any] = {}
        self.added: list[Any] = []
        self.commits: int = 0
        self.refreshed: list[Any] = []
        self.last_filters: list[Any] = []

    def add(self, obj: Any):
        self.added.append(obj)
        # Assign id if needed
        if getattr(obj, "id", None) is None:
            obj.id = uuid4()
        self.notifications.append(obj)
        self.by_id[obj.id] = obj

    def commit(self):
        self.commits += 1

    def refresh(self, obj: Any):
        self.refreshed.append(obj)

    def query(self, _model: Any) -> _Query:
        return _Query(self)


def test_create_notification_happy(monkeypatch):
    db = _DB()
    user_id = uuid4()

    # Patch ensure_uuid to produce UUID consistently
    monkeypatch.setattr(ns, "ensure_uuid", lambda s: UUID(str(s)))

    # Patch Notification constructor used by the service
    @dataclass
    class N:
        user_id: UUID
        title: str
        message: Optional[str]
        notification_type: str
        is_read: bool
        id: Optional[UUID] = None

    monkeypatch.setattr(ns, "Notification", N)

    created = ns.create_notification(str(user_id), "T", message="M", notification_type="info", db=db)
    assert created.title == "T"
    assert created.is_read is False
    assert db.commits == 1
    assert db.refreshed == [created]


def test_get_unread_notifications_and_count(monkeypatch):
    db = _DB()
    now = datetime.utcnow()

    # Provide real objects in db list (we won't use SQLAlchemy filtering)
    class RealN:
        def __init__(self, id: UUID, is_read: bool, created_at: datetime):
            self.id = id
            self.is_read = is_read
            self.created_at = created_at

    # Seed: two unread, one read
    db.notifications = [
        RealN(uuid4(), False, now - timedelta(seconds=10)),
        RealN(uuid4(), True, now),
        RealN(uuid4(), False, now),
    ]

    # Patch Notification columns so query-building doesn't crash
    class Model:
        user_id = _Col("user_id")
        is_read = _Col("is_read")
        created_at = _Col("created_at")

    monkeypatch.setattr(ns, "Notification", Model)
    monkeypatch.setattr(ns, "ensure_uuid", lambda s: UUID(str(uuid4())))  # value irrelevant to our FakeQuery

    unread = ns.get_unread_notifications(str(uuid4()), db, limit=10)
    assert [n.created_at for n in unread] == sorted([n.created_at for n in unread], reverse=True)
    assert len(unread) == 2
    assert ns.get_unread_count(str(uuid4()), db) == 2


def test_mark_as_read_found_and_not_found(monkeypatch):
    db = _DB()
    notif_id = uuid4()

    class RealN:
        def __init__(self, id: UUID, is_read: bool):
            self.id = id
            self.is_read = is_read

    n = RealN(notif_id, False)
    db.notifications = [n]
    db.by_id = {notif_id: n}

    class Model:
        id = _Col("id")

    monkeypatch.setattr(ns, "Notification", Model)
    monkeypatch.setattr(ns, "ensure_uuid", lambda s: UUID(str(s)))

    ns.mark_as_read(str(notif_id), db)
    assert n.is_read is True
    assert db.commits == 1

    # Not found should not commit
    before = db.commits
    ns.mark_as_read(str(uuid4()), db)
    assert db.commits == before


def test_mark_all_read_updates_only_unread(monkeypatch):
    db = _DB()

    class RealN:
        def __init__(self, id: UUID, is_read: bool):
            self.id = id
            self.is_read = is_read

    db.notifications = [RealN(uuid4(), False), RealN(uuid4(), True), RealN(uuid4(), False)]

    class Model:
        user_id = _Col("user_id")
        is_read = _Col("is_read")

    monkeypatch.setattr(ns, "Notification", Model)
    monkeypatch.setattr(ns, "ensure_uuid", lambda s: UUID(str(uuid4())))

    ns.mark_all_read(str(uuid4()), db)
    assert all(n.is_read for n in db.notifications)
    assert db.commits == 1
    # Ensure we used the correct SQL-ish filter, not a literal False
    assert ("is_read", "is", False) in db.last_filters


