from __future__ import annotations

import uuid

from amprenta_rag.automation import actions as act


class FakeModel:
    def __init__(self, **fields):
        self.__dict__.update(fields)
        self.id = fields.get("id", uuid.uuid4())


class FakeDB:
    def __init__(self):
        self.added = []
        self.refreshed = []

    def add(self, obj):
        self.added.append(obj)

    def commit(self):
        return None

    def refresh(self, obj):
        self.refreshed.append(obj)


class FakeIssue:
    def __init__(self, field: str, issue: str, severity: str):
        self.field = field
        self.issue = issue
        self.severity = severity


def test_send_notification_happy_path(monkeypatch) -> None:
    db = FakeDB()
    fake_note = FakeModel(id=uuid.uuid4(), user_id=uuid.uuid4())

    monkeypatch.setattr(act, "Notification", lambda **kwargs: fake_note)

    result = act.send_notification(
        {"user_id": str(uuid.uuid4()), "title": "Hello {name}", "message": "Msg {name}"},
        {"name": "World"},
        db,
    )

    assert result["status"] == "sent"
    assert db.added


def test_send_notification_requires_user(monkeypatch) -> None:
    db = FakeDB()
    try:
        act.send_notification({}, {}, db)
    except ValueError as e:
        assert "user_id" in str(e)


def test_add_note_and_errors(monkeypatch) -> None:
    db = FakeDB()
    fake_note = FakeModel(id=uuid.uuid4())
    monkeypatch.setattr(act, "Note", lambda **kwargs: fake_note)

    res = act.add_note(
        {"entity_type": "experiment", "entity_id": str(uuid.uuid4()), "content": "Hi {x}"},
        {"x": "there", "created_by_id": str(uuid.uuid4())},
        db,
    )
    assert res["status"] == "created"

    try:
        act.add_note({"entity_type": "exp"}, {}, db)
    except ValueError as e:
        assert "entity_type and entity_id" in str(e)


def test_run_validation(monkeypatch) -> None:
    issues = [FakeIssue("f1", "bad", "error"), FakeIssue("f2", "warn", "warning")]
    monkeypatch.setattr(act, "validate_experiment", lambda ent_id, db: issues)
    monkeypatch.setattr(act, "validate_compound", lambda ent_id, db: issues[:1])

    res = act.run_validation({"entity_type": "experiment", "entity_id": "123"}, {}, None)
    assert res["error_count"] == 1
    assert res["warning_count"] == 1

    res2 = act.run_validation({"entity_type": "compound", "entity_id": "123"}, {}, None)
    assert res2["error_count"] == 1

    try:
        act.run_validation({"entity_type": "other", "entity_id": "123"}, {}, None)
    except ValueError as e:
        assert "Unsupported" in str(e)

