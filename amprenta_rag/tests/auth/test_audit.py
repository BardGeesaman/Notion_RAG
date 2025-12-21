from __future__ import annotations
import pytest
from uuid import uuid4
from amprenta_rag.auth import audit

class FakeDB:
    def __init__(self):
        self.added = []
        
    def add(self, obj):
        self.added.append(obj)
        
    def commit(self):
        pass

def test_log_action_success(monkeypatch):
    fake_db = FakeDB()
    monkeypatch.setattr(audit, "get_db", lambda: iter([fake_db]))
    
    uid = str(uuid4())
    audit.log_action("test_action", user_id=uid, username="user", entity_type="exp", entity_id="123")
    
    assert len(fake_db.added) == 1
    entry = fake_db.added[0]
    assert entry.action == "test_action"
    assert entry.username == "user"
    assert entry.entity_type == "exp"

def test_log_action_error_handling(monkeypatch, caplog):
    # Simulate DB error
    def crash():
        raise RuntimeError("DB failed")
        yield
    
    monkeypatch.setattr(audit, "get_db", crash)
    
    audit.log_action("action")
    assert "Failed to log action" in caplog.text

def test_log_wrappers(monkeypatch):
    calls = []
    def fake_log_action(action, user_id=None, username=None, entity_type=None, entity_id=None, details=None, ip_address=None):
        calls.append(action)
    
    monkeypatch.setattr(audit, "log_action", fake_log_action)
    
    audit.log_login("u1", "user")
    audit.log_logout("u1", "user")
    audit.log_create("u1", "user", "type", "id")
    audit.log_update("u1", "user", "type", "id")
    audit.log_delete("u1", "user", "type", "id")
    
    assert calls == ["login", "logout", "create", "update", "delete"]

