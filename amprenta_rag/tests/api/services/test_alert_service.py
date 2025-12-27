from __future__ import annotations

from datetime import datetime, timedelta
from uuid import uuid4

from amprenta_rag.api.services import alert_service


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def desc(self):
        return ("desc", self.name)

    def is_(self, val):
        return ("is_val", self.name, val)

    def __hash__(self):
        return hash(self.name)


class FakeRepositoryNotification:
    id = Field("id")
    created_at = Field("created_at")
    is_read = Field("is_read")
    subscription_id = None
    dataset_id = None

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.user_id = kwargs.get("user_id")
        self.title = kwargs.get("title")
        self.description = kwargs.get("description")
        self.severity = kwargs.get("severity", "info")
        self.created_at = kwargs.get("created_at", datetime.utcnow())
        self.is_read = kwargs.get("is_read", False)
        self.subscription_id = kwargs.get("subscription_id")
        self.dataset_id = kwargs.get("dataset_id")


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, predicate):
        filtered = []
        for obj in self._data:
            ok = True
            if isinstance(predicate, tuple):
                kind = predicate[0]
                if kind == "eq":
                    field, val = predicate[1], predicate[2]
                    if getattr(obj, field, None) != val:
                        ok = False
                elif kind == "is_val":
                    field, val = predicate[1], predicate[2]
                    if getattr(obj, field, None) is not val:
                        ok = False
            if predicate is False:
                ok = False
            if ok:
                filtered.append(obj)
        return FakeQuery(filtered)

    def join(self, *args, **kwargs):
        return self

    def update(self, *args, **kwargs):
        return len(self._data)

    def order_by(self, *args, **kwargs):
        return self

    def limit(self, n):
        return FakeQuery(self._data[:n])

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self, alerts=None):
        self.alerts = alerts or []
        self.deleted = []

    def query(self, model):
        return FakeQuery(self.alerts)

    def commit(self):
        return None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_list_alerts_limits_and_sorts(monkeypatch):
    old = FakeRepositoryNotification(created_at=datetime.utcnow() - timedelta(days=1))
    new = FakeRepositoryNotification(created_at=datetime.utcnow())
    db = FakeSession([new, old])
    monkeypatch.setattr(alert_service, "db_session", lambda: db)
    monkeypatch.setattr(alert_service, "RepositoryNotification", FakeRepositoryNotification)

    rows = alert_service.list_alerts(user_id=None, unread_only=False)
    assert rows[0]["id"] == new.id


def test_mark_read(monkeypatch):
    alert = FakeRepositoryNotification(is_read=False)
    db = FakeSession([alert])
    monkeypatch.setattr(alert_service, "db_session", lambda: db)
    monkeypatch.setattr(alert_service, "RepositoryNotification", FakeRepositoryNotification)
    assert alert_service.mark_read(alert.id) is True
    assert alert.is_read is True
    assert alert_service.mark_read(uuid4()) is False


def test_delete_alert(monkeypatch):
    alert = FakeRepositoryNotification(is_read=False)
    db = FakeSession([alert])
    monkeypatch.setattr(alert_service, "db_session", lambda: db)
    monkeypatch.setattr(alert_service, "RepositoryNotification", FakeRepositoryNotification)
    updated = alert_service.mark_all_read(user_id=None)
    assert updated == 1

