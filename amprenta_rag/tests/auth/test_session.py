from __future__ import annotations

import time
from amprenta_rag.auth import session

def test_get_current_user(monkeypatch):
    monkeypatch.setattr(session.st, "session_state", {"user": {"name": "test"}})
    assert session.get_current_user() == {"name": "test"}
    
    monkeypatch.setattr(session.st, "session_state", {})
    assert session.get_current_user() is None

def test_set_current_user(monkeypatch):
    state = {}
    monkeypatch.setattr(session.st, "session_state", state)
    session.set_current_user({"name": "new"})
    assert state["user"] == {"name": "new"}

def test_clear_session(monkeypatch):
    state = {"user": "exists"}
    monkeypatch.setattr(session.st, "session_state", state)
    session.clear_session()
    assert "user" not in state

def test_is_authenticated(monkeypatch):
    monkeypatch.setattr(session.st, "session_state", {"user": "exists"})
    assert session.is_authenticated() is True
    
    monkeypatch.setattr(session.st, "session_state", {})
    assert session.is_authenticated() is False

def test_require_auth(monkeypatch):
    # Mirrors is_authenticated for now
    monkeypatch.setattr(session.st, "session_state", {"user": "exists"})
    assert session.require_auth() is True

def test_update_last_activity(monkeypatch):
    state = {}
    monkeypatch.setattr(session.st, "session_state", state)
    session.update_last_activity()
    assert "last_activity" in state
    assert isinstance(state["last_activity"], float)

def test_check_session_timeout(monkeypatch):
    state = {"user": "u"}
    monkeypatch.setattr(session.st, "session_state", state)
    
    # 1. Not authenticated
    monkeypatch.setattr(session.st, "session_state", {})
    assert session.check_session_timeout() is False
    
    # 2. Authenticated, no last activity -> initializes it, not timed out
    state = {"user": "u"}
    monkeypatch.setattr(session.st, "session_state", state)
    assert session.check_session_timeout() is False
    assert "last_activity" in state
    
    # 3. Authenticated, recent activity -> not timed out
    state["last_activity"] = time.time()
    assert session.check_session_timeout(timeout_minutes=30) is False
    
    # 4. Authenticated, old activity -> timed out
    state["last_activity"] = time.time() - (31 * 60) # 31 mins ago
    assert session.check_session_timeout(timeout_minutes=30) is True

def test_get_session_remaining(monkeypatch):
    # 1. Not authenticated
    monkeypatch.setattr(session.st, "session_state", {})
    assert session.get_session_remaining() == 0
    
    # 2. Authenticated, no activity
    state = {"user": "u"}
    monkeypatch.setattr(session.st, "session_state", state)
    assert session.get_session_remaining(30) == 30
    
    # 3. Authenticated, 10 mins elapsed
    state["last_activity"] = time.time() - (10 * 60)
    # roughly 20 mins remaining
    rem = session.get_session_remaining(30)
    assert 19 <= rem <= 20
    
    # 4. Authenticated, timed out
    state["last_activity"] = time.time() - (40 * 60)
    assert session.get_session_remaining(30) == 0
