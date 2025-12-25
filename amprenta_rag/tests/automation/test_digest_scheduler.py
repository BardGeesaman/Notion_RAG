from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import UUID, uuid4

import pytest

# Graceful skipping in minimal test environments
papermill = pytest.importorskip("papermill")
apscheduler = pytest.importorskip("apscheduler")
sklearn = pytest.importorskip("sklearn")

from amprenta_rag.automation import digest_scheduler as ds
from amprenta_rag.automation import notifications as notif


class _Query:
    def __init__(self, first_obj: Any = None):
        self._first = first_obj

    def filter(self, *_args: Any, **_kwargs: Any) -> "_Query":
        return self

    def first(self) -> Any:
        return self._first


class _DB:
    def __init__(self):
        self.added: list[Any] = []
        self.commits: int = 0
        self._schedule: Any = None
        self._program: Any = None

    def add(self, obj: Any) -> None:
        self.added.append(obj)
        if obj.__class__.__name__ == "DigestSchedule":
            self._schedule = obj

    def commit(self) -> None:
        self.commits += 1

    def refresh(self, _obj: Any) -> None:
        return None

    def query(self, model: Any) -> _Query:
        name = getattr(model, "__name__", str(model))
        if name.endswith("DigestSchedule"):
            return _Query(first_obj=self._schedule)
        if name.endswith("Program"):
            return _Query(first_obj=self._program)
        return _Query(first_obj=None)


class _DBCtx:
    def __init__(self, db: _DB):
        self._db = db

    def __enter__(self) -> _DB:
        return self._db

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False


@dataclass
class _Program:
    id: UUID
    name: str = "MyProgram"


def test_schedule_weekly_digest(monkeypatch, tmp_path):
    # Avoid registry/FS checks
    monkeypatch.setattr(ds, "_validate_notebook_path", lambda p: tmp_path / "x.ipynb")

    db = _DB()
    db._program = _Program(id=uuid4(), name="P")
    monkeypatch.setattr(ds, "db_session", lambda: _DBCtx(db))

    called: list[UUID] = []
    monkeypatch.setattr(ds.DigestScheduler, "add_digest_schedule", lambda self, sid: called.append(sid))

    sched = ds.DigestScheduler().schedule_weekly_digest(
        program_id=db._program.id,
        notebook_path="deploy/jupyterhub/templates/dose_response.ipynb",
        recipients=["a@x.com"],
        day_of_week="mon",
        hour=9,
    )

    assert sched.program_id == db._program.id
    assert sched.notebook_path.endswith(".ipynb")
    assert sched.schedule_cron == "0 9 * * mon"
    assert called == [UUID(str(sched.id))]


def test_execute_digest_papermill(monkeypatch, tmp_path):
    schedule_id = uuid4()
    program_id = uuid4()

    in_nb = tmp_path / "in.ipynb"
    in_nb.write_text("{}", encoding="utf-8")

    monkeypatch.setattr(ds, "_validate_notebook_path", lambda p: in_nb)

    calls: list[dict[str, Any]] = []

    class _PM:
        @staticmethod
        def execute_notebook(inp: str, outp: str, parameters: dict) -> None:
            calls.append({"inp": inp, "outp": outp, "parameters": parameters})
            Path(outp).write_text("{}", encoding="utf-8")

    class _NBFormat:
        @staticmethod
        def read(_p: str, as_version: int = 4) -> dict:
            return {"cells": [], "metadata": {}, "nbformat": as_version, "nbformat_minor": 0}

    class _Exporter:
        def from_notebook_node(self, _nb: dict):
            return "<html>ok</html>", {}

    monkeypatch.setattr(ds, "pm", _PM)
    monkeypatch.setattr(ds, "PAPERMILL_AVAILABLE", True)
    monkeypatch.setattr(ds, "nbformat", _NBFormat)
    monkeypatch.setattr(ds, "HTMLExporter", _Exporter)
    monkeypatch.setattr(ds, "NBCONVERT_AVAILABLE", True)

    res = ds.DigestScheduler()._execute(schedule_id=schedule_id, program_id=program_id, notebook_path="x.ipynb")
    assert res.status == "success"
    assert res.output_ipynb.exists()
    assert res.output_html and res.output_html.exists()
    assert calls and calls[0]["parameters"]["program_id"] == str(program_id)


def test_send_notification(monkeypatch):
    emailed: list[tuple[list[str], str, str]] = []
    slacked: list[str] = []

    monkeypatch.setattr(notif, "_send_email", lambda to, subject, body: emailed.append((to, subject, body)))
    monkeypatch.setattr(notif, "_send_slack", lambda message, webhook_url: slacked.append(message))

    out = notif.send_digest_notification(
        recipients=["a@x.com", "not-an-email"],
        digest_url="http://example/digest.html",
        program_name="Prog",
        slack_webhook_url="http://slack-webhook",
    )

    assert out["emailed"] == ["a@x.com"]
    assert emailed and "Weekly Prog digest is ready" in emailed[0][1]
    assert slacked and "Weekly Prog digest is ready" in slacked[0]


