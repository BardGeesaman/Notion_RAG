from __future__ import annotations

import pytest
from uuid import uuid4
from unittest.mock import MagicMock
from amprenta_rag.auth import audit

def test_log_action_success(monkeypatch):
    db_mock = MagicMock()
    
    # Mock generator for get_db
    def fake_get_db():
        yield db_mock
        
    monkeypatch.setattr(audit, "get_db", fake_get_db)
    
    audit.log_action(
        action="test",
        user_id=str(uuid4()),
        username="user",
        entity_type="type",
        entity_id=str(uuid4()),
        details={"foo": "bar"}
    )
    
    db_mock.add.assert_called_once()
    db_mock.commit.assert_called_once()
    
    # Verify log entry content
    entry = db_mock.add.call_args[0][0]
    assert entry.action == "test"
    assert entry.username == "user"
    assert entry.details == {"foo": "bar"}

def test_log_action_error_handled(monkeypatch, caplog):
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
