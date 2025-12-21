from __future__ import annotations

from amprenta_rag.eval import metrics


def test_evaluation_metrics_pass_rate_zero_and_nonzero():
    m0 = metrics.EvaluationMetrics(total_tasks=0, passed=0, failed=0, cost_tokens=0, p95_latency_ms=0.0, violations=0)
    assert m0.pass_rate == 0.0

    m1 = metrics.EvaluationMetrics(total_tasks=4, passed=3, failed=1, cost_tokens=10, p95_latency_ms=1.0, violations=0)
    assert m1.pass_rate == 0.75


def test_aggregate_metrics_computes_counts_and_p95():
    rows = [
        {"passed": 1, "latency_ms": 10, "cost_tokens": 5, "violations": 1},
        {"passed": 0, "latency_ms": 20, "cost_tokens": 3, "violations": 0},
        {"passed": 1, "latency_ms": 30, "cost_tokens": 2, "violations": 0},
    ]
    agg = metrics.aggregate_metrics(rows)
    assert agg.total_tasks == 3
    assert agg.passed == 2
    assert agg.failed == 1
    assert agg.cost_tokens == 10
    assert agg.violations == 1
    assert agg.p95_latency_ms == 30.0  # sorted: [10,20,30] idx=2

