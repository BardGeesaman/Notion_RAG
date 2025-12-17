from __future__ import annotations

import importlib
import importlib.util
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4, UUID
from unittest.mock import Mock

import pytest


class _Col:
    """Tiny SQLAlchemy-ish column stub to support == and .desc()."""

    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other):  # noqa: DunderEq
        return ("eq", self.name, other)

    def desc(self):
        return ("desc", self.name)


class FakeReportArtifact:
    id = _Col("id")
    entity_type = _Col("entity_type")
    entity_id = _Col("entity_id")
    params_hash = _Col("params_hash")
    created_at = _Col("created_at")

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        # mimic SQLAlchemy id default; tests can override
        if not hasattr(self, "id"):
            self.id = uuid4()


class FakeQuery:
    def __init__(self, db, first_queue=None, all_result=None):
        self._db = db
        self._first_queue = list(first_queue or [])
        self._all_result = all_result
        self.filter_calls = []
        self.order_by_calls = []
        self.limit_calls = []

    def filter(self, *args, **kwargs):
        self.filter_calls.append((args, kwargs))
        return self

    def order_by(self, *args, **kwargs):
        self.order_by_calls.append((args, kwargs))
        return self

    def limit(self, n: int):
        self.limit_calls.append(n)
        return self

    def first(self):
        if self._first_queue:
            return self._first_queue.pop(0)
        return None

    def all(self):
        return list(self._all_result or [])


class FakeDB:
    def __init__(self, *, first_queue=None, all_result=None):
        self._first_queue = list(first_queue or [])
        self._all_result = list(all_result or [])
        self.add = Mock()
        self.delete = Mock()
        self.commit = Mock()
        self.refresh = Mock()

    def query(self, _model):
        return FakeQuery(self, first_queue=self._first_queue, all_result=self._all_result)


@pytest.mark.unit
def test_artifact_registry_save_get_list_cached_delete(monkeypatch):
    from amprenta_rag.reports import artifact_registry as ar

    # Patch model class to avoid SQLAlchemy dependency surface.
    monkeypatch.setattr(ar, "ReportArtifact", FakeReportArtifact)

    entity_id = uuid4()
    user_id = uuid4()
    created = FakeReportArtifact(
        entity_type="program",
        entity_id=entity_id,
        format="html",
        file_path="/tmp/r.html",
        params_hash="abc",
        created_by_id=user_id,
    )
    created.id = uuid4()

    # DB will be reused by multiple functions; each opens its own session.
    db_for_save = FakeDB()

    @contextmanager
    def fake_db_session_save():
        yield db_for_save

    monkeypatch.setattr(ar, "db_session", fake_db_session_save)

    out = ar.save_artifact(
        entity_type="program",
        entity_id=entity_id,
        format="html",
        file_path="/tmp/r.html",
        params_hash="abc",
        user_id=user_id,
    )
    assert isinstance(out, FakeReportArtifact)
    db_for_save.add.assert_called_once()
    db_for_save.commit.assert_called_once()
    db_for_save.refresh.assert_called_once_with(out)

    # get_artifact
    db_for_get = FakeDB(first_queue=[created])

    @contextmanager
    def fake_db_session_get():
        yield db_for_get

    monkeypatch.setattr(ar, "db_session", fake_db_session_get)
    got = ar.get_artifact(created.id)
    assert got is created

    # list_artifacts
    listed = [created]
    db_for_list = FakeDB(all_result=listed)

    @contextmanager
    def fake_db_session_list():
        yield db_for_list

    monkeypatch.setattr(ar, "db_session", fake_db_session_list)
    out_list = ar.list_artifacts("program", entity_id, limit=12)
    assert out_list == listed

    # get_cached_artifact
    db_for_cached = FakeDB(first_queue=[created])

    @contextmanager
    def fake_db_session_cached():
        yield db_for_cached

    monkeypatch.setattr(ar, "db_session", fake_db_session_cached)
    cached = ar.get_cached_artifact("program", entity_id, params_hash="abc")
    assert cached is created

    # delete_artifact true
    db_for_del = FakeDB(first_queue=[created])

    @contextmanager
    def fake_db_session_del():
        yield db_for_del

    monkeypatch.setattr(ar, "db_session", fake_db_session_del)
    assert ar.delete_artifact(created.id) is True
    db_for_del.delete.assert_called_once_with(created)
    db_for_del.commit.assert_called_once()

    # delete_artifact false
    db_for_del2 = FakeDB(first_queue=[None])

    @contextmanager
    def fake_db_session_del2():
        yield db_for_del2

    monkeypatch.setattr(ar, "db_session", fake_db_session_del2)
    assert ar.delete_artifact(uuid4()) is False
    db_for_del2.delete.assert_not_called()


@pytest.mark.unit
def test_scheduler_imports_and_functions_callable(monkeypatch):
    """
    Scheduler should import and core functions should be callable with dependencies mocked.
    """
    # Allow this test suite to run even if APScheduler isn't installed in the test env.
    # The runtime image may include it; here we only need to validate module structure.
    if "apscheduler" not in sys.modules:
        import types

        apscheduler_mod = types.ModuleType("apscheduler")
        triggers_mod = types.ModuleType("apscheduler.triggers")
        cron_mod = types.ModuleType("apscheduler.triggers.cron")

        class _CronTrigger:
            @staticmethod
            def from_crontab(expr: str):
                return ("cron", expr)

        cron_mod.CronTrigger = _CronTrigger

        monkeypatch.setitem(sys.modules, "apscheduler", apscheduler_mod)
        monkeypatch.setitem(sys.modules, "apscheduler.triggers", triggers_mod)
        monkeypatch.setitem(sys.modules, "apscheduler.triggers.cron", cron_mod)

    sched = importlib.import_module("amprenta_rag.reports.scheduler")

    assert callable(getattr(sched, "run_scheduled_report"))
    assert callable(getattr(sched, "add_report_schedule"))
    assert callable(getattr(sched, "remove_report_schedule"))
    assert callable(getattr(sched, "load_report_schedules"))


@pytest.mark.unit
def test_scheduler_has_no_duplicate_definitions_in_source():
    """
    Guardrail: scheduler.py should not contain duplicate function definitions.
    This currently fails if the module body was accidentally duplicated.
    """
    path = Path(__file__).resolve().parents[3] / "amprenta_rag/reports/scheduler.py"
    text = path.read_text(encoding="utf-8")
    assert text.count("def run_scheduled_report") == 1
    assert text.count("def add_report_schedule") == 1
    assert text.count("def remove_report_schedule") == 1
    assert text.count("def load_report_schedules") == 1


def _load_module_from_path(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


@pytest.mark.unit
def test_report_history_uuid_validation_rejects_invalid(monkeypatch):
    """
    report_history._render_manual_generate() should reject invalid UUID input.
    We stub streamlit and drive the form submission path.
    """
    errors = []

    class _Form:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class FakeStreamlit:
        def subheader(self, *a, **k):
            return None

        def form(self, key):
            return _Form()

        def selectbox(self, *a, **k):
            # entity_type then fmt
            if "Entity type" in a[0]:
                return "program"
            return "html"

        def text_input(self, *a, **k):
            return "not-a-uuid"

        def form_submit_button(self, *a, **k):
            return True

        def error(self, msg):
            errors.append(msg)

        # used elsewhere but not in this path
        def info(self, *a, **k):
            return None

        def title(self, *a, **k):
            return None

        def markdown(self, *a, **k):
            return None

    # Ensure our fake streamlit is used during module import
    monkeypatch.setitem(sys.modules, "streamlit", FakeStreamlit())
    # Stub scheduler import to avoid failures if scheduler.py is broken; UUID validation should be testable in isolation.
    monkeypatch.setitem(
        sys.modules,
        "amprenta_rag.reports.scheduler",
        SimpleNamespace(add_report_schedule=lambda *a, **k: None, remove_report_schedule=lambda *a, **k: None),
    )

    mod = _load_module_from_path(
        "report_history",
        Path(__file__).resolve().parents[3] / "scripts/dashboard/pages/report_history.py",
    )

    mod._render_manual_generate()
    assert errors
    assert "Invalid entity ID format" in errors[-1]


