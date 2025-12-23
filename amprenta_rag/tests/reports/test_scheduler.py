from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional
from uuid import UUID, uuid4


from amprenta_rag.reports import scheduler as rs


class _ReportScheduleModel:
    # Only used to satisfy attribute access in filter expressions.
    id = object()
    enabled = True


class _Query:
    def __init__(self, first_obj: Any = None, all_objs: Optional[List[Any]] = None):
        self._first = first_obj
        self._all = list(all_objs or [])

    def filter(self, *_args: Any, **_kwargs: Any) -> "_Query":
        return self

    def first(self) -> Any:
        return self._first

    def all(self) -> List[Any]:
        return list(self._all)


class _DB:
    def __init__(self, first_obj: Any = None, all_objs: Optional[List[Any]] = None):
        self._q = _Query(first_obj=first_obj, all_objs=all_objs)
        self.commits: int = 0

    def query(self, _model: Any) -> _Query:
        return self._q

    def commit(self) -> None:
        self.commits += 1


class _DBCtx:
    def __init__(self, db: _DB):
        self._db = db

    def __enter__(self) -> _DB:
        return self._db

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False


@dataclass
class _Schedule:
    id: UUID
    name: str = "Sched"
    enabled: bool = True
    entity_type: str = "dataset"
    entity_id: Optional[str] = None
    format: str = "pdf"
    created_by_id: Optional[str] = None
    cron_expression: str = "0 0 * * *"
    last_run_at: Any = None


class _FakeScheduler:
    def __init__(self):
        self.jobs: dict[str, dict[str, Any]] = {}
        self.removed: list[str] = []

    def get_job(self, job_id: str):
        return self.jobs.get(job_id)

    def remove_job(self, job_id: str) -> None:
        self.removed.append(job_id)
        self.jobs.pop(job_id, None)

    def add_job(self, func, trigger, args, id: str, name: str, replace_existing: bool):
        self.jobs[id] = {
            "func": func,
            "trigger": trigger,
            "args": args,
            "name": name,
            "replace_existing": replace_existing,
        }


def test_params_hash_stable():
    eid = uuid4()
    h1 = rs._params_hash("dataset", eid, "pdf")
    h2 = rs._params_hash("dataset", eid, "pdf")
    h3 = rs._params_hash("dataset", eid, "html")
    assert h1 == h2
    assert h1 != h3


def test_run_scheduled_report_not_found(monkeypatch):
    db = _DB(first_obj=None)
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)

    called: list[str] = []
    monkeypatch.setattr(rs, "generate_report", lambda *a, **k: called.append("gen"))

    rs.run_scheduled_report(uuid4())
    assert called == []


def test_run_scheduled_report_disabled(monkeypatch):
    sched = _Schedule(id=uuid4(), enabled=False)
    db = _DB(first_obj=sched)
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)

    called: list[str] = []
    monkeypatch.setattr(rs, "generate_report", lambda *a, **k: called.append("gen"))

    rs.run_scheduled_report(sched.id)
    assert called == []


def test_run_scheduled_report_success_saves_artifact(monkeypatch, tmp_path):
    eid = uuid4()
    sched = _Schedule(
        id=uuid4(),
        enabled=True,
        entity_type="dataset",
        entity_id=str(eid),
        format="pdf",
        created_by_id=str(uuid4()),
    )
    db = _DB(first_obj=sched)
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)
    monkeypatch.setattr(rs, "ensure_uuid", lambda x: UUID(str(x)) if x else None)

    out_path = tmp_path / "r.pdf"
    out_path.write_bytes(b"x")
    monkeypatch.setattr(rs, "generate_report", lambda *a, **k: out_path)

    saved: list[dict[str, Any]] = []

    def fake_save_artifact(**kwargs):
        saved.append(kwargs)

    monkeypatch.setattr(rs, "save_artifact", fake_save_artifact)

    rs.run_scheduled_report(sched.id)

    assert db.commits == 1
    assert sched.last_run_at is not None
    assert len(saved) == 1
    assert saved[0]["entity_type"] == "dataset"
    assert saved[0]["entity_id"] == eid
    assert saved[0]["format"] == "pdf"
    assert saved[0]["file_path"] == out_path
    assert saved[0]["user_id"] is not None
    assert saved[0]["params_hash"]


def test_run_scheduled_report_generation_error(monkeypatch):
    eid = uuid4()
    sched = _Schedule(id=uuid4(), enabled=True, entity_id=str(eid), created_by_id=None)
    db = _DB(first_obj=sched)
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)
    monkeypatch.setattr(rs, "ensure_uuid", lambda x: UUID(str(x)) if x else None)

    monkeypatch.setattr(rs, "generate_report", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))

    saved: list[bool] = []
    monkeypatch.setattr(rs, "save_artifact", lambda **k: saved.append(True))

    rs.run_scheduled_report(sched.id)
    assert saved == []


def test_add_report_schedule_adds_job(monkeypatch):
    sched_id = uuid4()
    sched = _Schedule(id=sched_id, enabled=True, name="MySched", cron_expression="*/5 * * * *")
    db = _DB(first_obj=sched)
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)

    fake_sched = _FakeScheduler()
    # Pretend an old job exists
    fake_sched.jobs[f"report_{sched_id}"] = {"old": True}
    monkeypatch.setattr(rs, "get_scheduler", lambda: fake_sched)
    monkeypatch.setattr(rs, "start_scheduler", lambda: None)

    class _Trig:
        pass

    called: list[str] = []
    monkeypatch.setattr(
        rs.CronTrigger,
        "from_crontab",
        lambda expr: called.append(expr) or _Trig(),
    )

    rs.add_report_schedule(sched_id)

    assert called == ["*/5 * * * *"]
    assert f"report_{sched_id}" in fake_sched.jobs
    assert fake_sched.removed == [f"report_{sched_id}"]


def test_remove_report_schedule(monkeypatch):
    sched_id = uuid4()
    fake_sched = _FakeScheduler()
    fake_sched.jobs[f"report_{sched_id}"] = {"old": True}
    monkeypatch.setattr(rs, "get_scheduler", lambda: fake_sched)

    rs.remove_report_schedule(sched_id)
    assert fake_sched.removed == [f"report_{sched_id}"]
    assert fake_sched.get_job(f"report_{sched_id}") is None


def test_load_report_schedules_calls_add(monkeypatch):
    s1 = _Schedule(id=uuid4(), enabled=True, cron_expression="0 1 * * *")
    s2 = _Schedule(id=uuid4(), enabled=True, cron_expression="0 2 * * *")
    db = _DB(all_objs=[s1, s2])
    monkeypatch.setattr(rs, "db_session", lambda: _DBCtx(db))
    monkeypatch.setattr(rs, "ReportSchedule", _ReportScheduleModel)
    monkeypatch.setattr(rs, "ensure_uuid", lambda x: x)

    added: list[tuple[UUID, str]] = []
    monkeypatch.setattr(rs, "add_report_schedule", lambda sid, cron=None: added.append((sid, cron or "")))

    rs.load_report_schedules()
    assert added == [(s1.id, s1.cron_expression), (s2.id, s2.cron_expression)]


