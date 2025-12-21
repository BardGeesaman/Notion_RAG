from __future__ import annotations

import logging
import pytest
from dataclasses import dataclass
from datetime import datetime
from uuid import uuid4
from amprenta_rag.auth import audit


@dataclass
class _FakeAuditLog:
    user_id: object = None
    username: str | None = None
    action: str | None = None
    entity_type: str | None = None
    entity_id: str | None = None
    details: dict | None = None
    ip_address: str | None = None
    timestamp: datetime | None = None

def test_log_action_success(monkeypatch):
    added: list[_FakeAuditLog] = []
    committed: list[bool] = []
    closed: list[bool] = []

    class _DB:
        def add(self, obj):
            added.append(obj)

        def commit(self):
            committed.append(True)

    def fake_get_db():
        try:
            yield _DB()
        finally:
            closed.append(True)

    # Avoid touching ORM / real UUID parsing in unit tests.
    monkeypatch.setattr(audit, "get_db", fake_get_db)
    monkeypatch.setattr(audit, "AuditLog", _FakeAuditLog)
    monkeypatch.setattr(audit, "ensure_uuid", lambda s: s)

    user_id = str(uuid4())
    entity_id = str(uuid4())
    audit.log_action(
        action="test",
        user_id=user_id,
        username="user",
        entity_type="type",
        entity_id=entity_id,
        details={"foo": "bar"},
        ip_address="1.2.3.4",
    )

    assert committed == [True]
    assert closed == [True]
    assert len(added) == 1
    entry = added[0]
    assert entry.action == "test"
    assert entry.username == "user"
    assert entry.entity_type == "type"
    assert entry.entity_id == entity_id
    assert entry.details == {"foo": "bar"}
    assert entry.ip_address == "1.2.3.4"
    assert entry.user_id == user_id
    assert isinstance(entry.timestamp, datetime)


def test_log_action_user_id_test_special_case(monkeypatch):
    added: list[_FakeAuditLog] = []

    class _DB:
        def add(self, obj):
            added.append(obj)

        def commit(self):
            return None

    def fake_get_db():
        yield _DB()

    monkeypatch.setattr(audit, "get_db", fake_get_db)
    monkeypatch.setattr(audit, "AuditLog", _FakeAuditLog)

    # If ensure_uuid were called we'd fail the test; "test" should bypass it.
    monkeypatch.setattr(audit, "ensure_uuid", lambda _s: (_ for _ in ()).throw(AssertionError("ensure_uuid called")))

    audit.log_action(action="login", user_id="test", username="u")
    assert len(added) == 1
    assert added[0].user_id is None

def test_log_action_error_handled(monkeypatch, caplog):
    caplog.set_level(logging.WARNING)
    monkeypatch.setattr(audit, "get_db", lambda: (_ for _ in range(1)).throw(RuntimeError("Boom")))
    audit.log_action("action")
    assert "Failed to log action" in caplog.text

def test_log_action_wrappers(monkeypatch):
    calls = []
    
    def fake_log_action(action, *args, **kwargs):
        calls.append(action)
        
    monkeypatch.setattr(audit, "log_action", fake_log_action)
    
    audit.log_login("u", "user")
    assert calls[-1] == "login"
    
    audit.log_logout("u", "user")
    assert calls[-1] == "logout"
    
    audit.log_create("u", "user", "type", "id")
    assert calls[-1] == "create"
    
    audit.log_update("u", "user", "type", "id")
    assert calls[-1] == "update"
    
    audit.log_delete("u", "user", "type", "id")
    assert calls[-1] == "delete"
