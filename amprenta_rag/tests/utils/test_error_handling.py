from __future__ import annotations

import time

import pytest

from amprenta_rag.utils import error_handling


def test_retry_with_backoff_succeeds_after_retries(monkeypatch):
    calls = {"count": 0}
    monkeypatch.setattr(time, "sleep", lambda *a, **k: None)

    @error_handling.retry_with_backoff(error_handling.RetryConfig(max_attempts=3, initial_delay=0.01))
    def flaky():
        calls["count"] += 1
        if calls["count"] < 3:
            raise ValueError("fail")
        return "ok"

    assert flaky() == "ok"
    assert calls["count"] == 3


def test_retry_with_backoff_raises_after_max(monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda *a, **k: None)

    @error_handling.retry_with_backoff(error_handling.RetryConfig(max_attempts=2, initial_delay=0))
    def always_fail():
        raise RuntimeError("nope")

    with pytest.raises(RuntimeError):
        always_fail()


def test_circuit_breaker_half_open_recovery(monkeypatch):
    cb = error_handling.CircuitBreaker(failure_threshold=1, recovery_timeout=0, expected_exception=ValueError)

    with pytest.raises(ValueError):
        cb.call(lambda: (_ for _ in ()).throw(ValueError("fail")))

    # Next call should move to half-open (timeout 0) and succeed
    assert cb.call(lambda: "ok") == "ok"
    assert cb.state == "closed"


def test_graceful_degradation_fallback(monkeypatch):
    logs = []
    monkeypatch.setattr(error_handling, "logger", type("L", (), {"warning": lambda *a, **k: logs.append(a), "error": lambda *a, **k: logs.append(a)})())

    @error_handling.graceful_degradation(fallback_value="fallback")
    def boom():
        raise RuntimeError("x")

    assert boom() == "fallback"
    assert logs, "Expected warning log"


def test_error_tracker_records_and_caps():
    tracker = error_handling.ErrorTracker(max_errors=2)
    tracker.record_error("op", ValueError("a"))
    tracker.record_error("op", RuntimeError("b"))
    tracker.record_error("op", RuntimeError("c"))
    summary = tracker.get_error_summary()
    assert summary["total_errors"] == 2
    assert summary["error_types"]["RuntimeError"] == 2

