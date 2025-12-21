from __future__ import annotations

import time

import pytest

from amprenta_rag.utils import health_check


class FakeResp:
    def __init__(self, status_code=200, data=None):
        self.status_code = status_code
        self._data = data or {}

    def json(self):
        return self._data


def test_health_check_configuration_only(monkeypatch):
    hc = health_check.HealthChecker()

    fake_cfg = type("Cfg", (), {"notion": type("N", (), {"base_url": "", "api_key": "", "version": ""})(), "pinecone": type("P", (), {"index_name": "idx"})(), "openai": type("O", (), {"api_key": "k"})()})
    monkeypatch.setattr(health_check, "get_config", lambda: fake_cfg)
    class FakeRequests:
        @staticmethod
        def get(*a, **k):
            return FakeResp(401)

    monkeypatch.setitem(health_check.__dict__, "requests", FakeRequests)
    monkeypatch.setattr(health_check, "get_pinecone_index", lambda: None, raising=False)
    monkeypatch.setattr(health_check, "openai", type("O", (), {"OpenAI": lambda api_key: type("C", (), {"models": type("M", (), {"list": lambda self: []})()})()})(), raising=False)

    # Skip ingestion checks
    hc._check_notion()
    hc._check_pinecone()
    hc._check_openai()
    hc._check_configuration()

    result = hc.check_all()
    assert result["components"]["total"] == len(hc.results)
    assert result["overall_status"] in {"healthy", "degraded", "unhealthy"}


def test_result_to_dict_contains_fields():
    hc = health_check.HealthChecker()
    res = health_check.HealthCheckResult("comp", "healthy", "ok", latency_ms=1.23, details={"k": "v"})
    data = hc._result_to_dict(res)
    assert data["component"] == "comp"
    assert data["status"] == "healthy"
    assert data["latency_ms"] == 1.23


def test_check_all_aggregates(monkeypatch):
    hc = health_check.HealthChecker()

    monkeypatch.setattr(hc, "_check_notion", lambda: hc.results.append(health_check.HealthCheckResult("n", "healthy", "ok")))
    monkeypatch.setattr(hc, "_check_pinecone", lambda: hc.results.append(health_check.HealthCheckResult("p", "degraded", "slow")))
    monkeypatch.setattr(hc, "_check_openai", lambda: hc.results.append(health_check.HealthCheckResult("o", "unhealthy", "down")))
    monkeypatch.setattr(hc, "_check_configuration", lambda: hc.results.append(health_check.HealthCheckResult("c", "healthy", "ok")))

    summary = hc.check_all()
    assert summary["components"]["healthy"] == 2
    assert summary["components"]["degraded"] == 1
    assert summary["components"]["unhealthy"] == 1
    assert summary["components"]["total"] == 4

