from __future__ import annotations

from datetime import datetime
from uuid import uuid4

from amprenta_rag.api.services import subscription_service


class Field:
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):
        return ("eq", self.name, other)

    def desc(self):
        return ("desc", self.name)


class FakeSubscription:
    id = Field("id")
    created_at = Field("created_at")
    user_id = Field("user_id")
    user_id_field = Field("user_id")
    name_field = Field("name")

    def __init__(self, **kwargs):
        self.id = kwargs.get("id", uuid4())
        self.user_id = kwargs.get("user_id")
        self.name = kwargs.get("name", "sub")
        self.repository_source = kwargs.get("repository_source", "src")
        self.query_params = kwargs.get("query_params", {})
        self.notify_email = kwargs.get("notify_email", False)
        self.notify_in_app = kwargs.get("notify_in_app", False)
        self.is_active = kwargs.get("is_active", True)
        self.repository = kwargs.get("repository")
        self.frequency = kwargs.get("frequency")
        self.created_at = kwargs.get("created_at", datetime.utcnow())
        self.updated_at = kwargs.get("updated_at")
        self.last_checked = kwargs.get("last_checked")
        self.created_at_field = self.created_at


class FakeQuery:
    def __init__(self, data):
        self._data = list(data)

    def filter(self, predicate):
        if isinstance(predicate, tuple):
            _, field, value = predicate
            filtered = [s for s in self._data if getattr(s, field, None) == value]
            return FakeQuery(filtered)
        return self

    def order_by(self, *_):
        return self

    def all(self):
        return list(self._data)

    def first(self):
        return self._data[0] if self._data else None


class FakeSession:
    def __init__(self, subs=None):
        self.subs = subs or []
        self.deleted = []

    def query(self, model):
        return FakeQuery(self.subs)

    def add(self, obj):
        self.subs.append(obj)

    def commit(self):
        return None

    def refresh(self, obj):
        return None

    def delete(self, obj):
        self.deleted.append(obj)
        self.subs = [s for s in self.subs if s is not obj]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_list_subscriptions(monkeypatch):
    s1 = FakeSubscription(user_id="u1")
    s2 = FakeSubscription(user_id="u2")
    monkeypatch.setattr(subscription_service, "db_session", lambda: FakeSession([s1, s2]))
    monkeypatch.setattr(subscription_service, "RepositorySubscription", FakeSubscription)
    subs = subscription_service.list_subscriptions(user_id=None)
    assert len(subs) == 2


def test_get_subscriptions_for_user(monkeypatch):
    s1 = FakeSubscription(user_id="u1")
    s2 = FakeSubscription(user_id="u2")
    monkeypatch.setattr(subscription_service, "db_session", lambda: FakeSession([s1, s2]))
    monkeypatch.setattr(subscription_service, "RepositorySubscription", FakeSubscription)
    subs = subscription_service.list_subscriptions(user_id="u1")
    assert subs and subs[0]["user_id"] == "u1"


def test_delete_subscription(monkeypatch):
    s1 = FakeSubscription(user_id="u1")
    db = FakeSession([s1])
    monkeypatch.setattr(subscription_service, "db_session", lambda: db)
    monkeypatch.setattr(subscription_service, "RepositorySubscription", FakeSubscription)
    assert subscription_service.delete_subscription(s1.id, user_id=None) is True
    assert s1 not in db.subs
    assert subscription_service.delete_subscription(uuid4(), user_id=None) is False

