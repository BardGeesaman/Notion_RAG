from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace
from uuid import uuid4
from unittest.mock import Mock

import pytest
from fastapi.testclient import TestClient


@pytest.mark.unit
def test_scan_all_datasets_returns_dataset_quality_reports(monkeypatch):
    from amprenta_rag.analysis import quality_watcher as qw

    id1 = uuid4()
    id2 = uuid4()
    ds1 = SimpleNamespace(id=id1, name="DS1")
    ds2 = SimpleNamespace(id=id2, name=None)

    db = SimpleNamespace(query=lambda model: SimpleNamespace(all=lambda: [ds1, ds2]))

    @contextmanager
    def fake_db_session():
        yield db

    monkeypatch.setattr(qw, "db_session", fake_db_session)
    monkeypatch.setattr(
        qw,
        "compute_quality_score",
        Mock(
            side_effect=[
                {"score": 90.0, "status": "high", "issues": [], "metrics": {"m": 1.0}},
                {"score": 10.0, "status": "low", "issues": ["x"], "metrics": {}},
            ]
        ),
    )

    reports = qw.scan_all_datasets()
    assert len(reports) == 2
    assert reports[0].dataset_id == id1
    assert reports[0].dataset_name == "DS1"
    assert reports[0].score == 90.0
    assert reports[0].status == "high"
    assert reports[0].issues == []
    assert reports[0].metrics == {"m": 1.0}

    # name fallback to UUID string
    assert reports[1].dataset_id == id2
    assert reports[1].dataset_name == str(id2)
    assert reports[1].score == 10.0
    assert reports[1].status == "low"
    assert reports[1].issues == ["x"]
    assert reports[1].metrics == {}


@pytest.mark.unit
def test_get_quality_summary_returns_counts_by_status(monkeypatch):
    from amprenta_rag.analysis import quality_watcher as qw

    monkeypatch.setattr(
        qw,
        "scan_all_datasets",
        Mock(
            return_value=[
                qw.DatasetQualityReport(uuid4(), "a", 90, "high", [], {}),
                qw.DatasetQualityReport(uuid4(), "b", 60, "medium", [], {}),
                qw.DatasetQualityReport(uuid4(), "c", 10, "low", [], {}),
                qw.DatasetQualityReport(uuid4(), "d", 55, "weird", [], {}),  # treated as low
            ]
        ),
    )
    summary = qw.get_quality_summary()
    assert summary == {"high": 1, "medium": 1, "low": 2, "total": 4}


@pytest.mark.unit
def test_get_low_quality_datasets_filters_correctly(monkeypatch):
    from amprenta_rag.analysis import quality_watcher as qw

    r1 = qw.DatasetQualityReport(uuid4(), "a", 49.9, "low", [], {})
    r2 = qw.DatasetQualityReport(uuid4(), "b", 50.0, "medium", [], {})
    r3 = qw.DatasetQualityReport(uuid4(), "c", 0.0, "low", [], {})
    monkeypatch.setattr(qw, "scan_all_datasets", Mock(return_value=[r1, r2, r3]))

    assert qw.get_low_quality_datasets(threshold=50.0) == [r1, r3]
    assert qw.get_low_quality_datasets(threshold=1.0) == [r3]


@pytest.fixture()
def api_client():
    from amprenta_rag.api.main import app

    return TestClient(app)


@pytest.mark.api
def test_get_quality_summary_endpoint(api_client, monkeypatch):
    from amprenta_rag.api.routers import quality as quality_router

    monkeypatch.setattr(
        quality_router,
        "get_quality_summary",
        Mock(return_value={"high": 2, "medium": 1, "low": 3, "total": 6}),
    )

    resp = api_client.get("/api/v1/quality/summary")
    assert resp.status_code == 200
    assert resp.json() == {"total": 6, "high": 2, "medium": 1, "low": 3}


@pytest.mark.api
def test_get_quality_datasets_endpoint(api_client, monkeypatch):
    from amprenta_rag.api.routers import quality as quality_router
    from amprenta_rag.analysis.quality_watcher import DatasetQualityReport as QReport

    id1 = uuid4()
    id2 = uuid4()
    monkeypatch.setattr(
        quality_router,
        "scan_all_datasets",
        Mock(
            return_value=[
                QReport(id1, "DS1", 88.8, "high", [], {"m": 1.0}),
                QReport(id2, "DS2", 12.3, "low", ["missing"], {}),
            ]
        ),
    )

    resp = api_client.get("/api/v1/quality/datasets")
    assert resp.status_code == 200
    body = resp.json()
    assert body == [
        {
            "dataset_id": str(id1),
            "dataset_name": "DS1",
            "score": 88.8,
            "status": "high",
            "issues": [],
            "metrics": {"m": 1.0},
        },
        {
            "dataset_id": str(id2),
            "dataset_name": "DS2",
            "score": 12.3,
            "status": "low",
            "issues": ["missing"],
            "metrics": {},
        },
    ]


@pytest.mark.api
def test_get_quality_alerts_endpoint_filters_by_threshold(api_client, monkeypatch):
    from amprenta_rag.api.routers import quality as quality_router
    from amprenta_rag.analysis.quality_watcher import DatasetQualityReport as QReport

    id1 = uuid4()
    monkeypatch.setattr(
        quality_router,
        "get_low_quality_datasets",
        Mock(return_value=[QReport(id1, "DS1", 10.0, "low", ["x"], {})]),
    )

    resp = api_client.get("/api/v1/quality/alerts?threshold=50")
    assert resp.status_code == 200
    assert resp.json()[0]["dataset_id"] == str(id1)
    quality_router.get_low_quality_datasets.assert_called_once_with(threshold=50.0)


@pytest.mark.api
def test_get_quality_alerts_threshold_validation(api_client):
    resp = api_client.get("/api/v1/quality/alerts?threshold=101")
    assert resp.status_code == 422


