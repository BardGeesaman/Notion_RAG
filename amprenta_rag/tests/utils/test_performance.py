from __future__ import annotations

import time

from amprenta_rag.utils import performance


def test_performance_timer_elapsed_and_logging(monkeypatch):
    logs = []
    monkeypatch.setattr(performance, "logger", type("L", (), {"info": lambda *a, **k: logs.append(a)})())

    with performance.PerformanceTimer("op", log_threshold=0) as t:
        time.sleep(0.01)
    assert t.elapsed_time() > 0
    assert logs, "Expected log when threshold is 0"


def test_timed_decorator_wraps_function():
    @performance.timed("fn", log_threshold=0)
    def add(a, b):
        return a + b

    result = add(2, 3)
    assert result == 5


def test_timer_context_manager():
    with performance.timer("ctx", log_threshold=0):
        x = 1 + 1
    assert x == 2


def test_performance_metrics_record_and_summary():
    metrics = performance.PerformanceMetrics()
    metrics.record("op", 1.0)
    metrics.record("op", 2.0)
    stats = metrics.get_statistics("op")
    assert stats and stats["count"] == 2 and stats["total"] == 3.0

    summary = metrics.get_summary()
    assert "op" in summary

    metrics.reset()
    assert metrics.get_statistics("op") is None

