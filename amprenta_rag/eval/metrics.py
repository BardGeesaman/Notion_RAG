from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass(frozen=True)
class EvaluationMetrics:
    total_tasks: int
    passed: int
    failed: int
    cost_tokens: int
    p95_latency_ms: float
    violations: int  # number of hard constraint violations

    @property
    def pass_rate(self) -> float:
        if self.total_tasks == 0:
            return 0.0
        return self.passed / self.total_tasks


def aggregate_metrics(rows: List[Dict[str, float | int]]) -> EvaluationMetrics:
    total = len(rows)
    passed = sum(1 for r in rows if r.get("passed") == 1)
    failed = total - passed
    cost = int(sum(int(r.get("cost_tokens", 0)) for r in rows))
    violations = int(sum(int(r.get("violations", 0)) for r in rows))

    # Robust p95 for small samples
    latencies = sorted(float(r.get("latency_ms", 0.0)) for r in rows)
    if latencies:
        idx = max(0, min(len(latencies) - 1, int(0.95 * (len(latencies) - 1))))
        p95 = latencies[idx]
    else:
        p95 = 0.0

    return EvaluationMetrics(
        total_tasks=total,
        passed=passed,
        failed=failed,
        cost_tokens=cost,
        p95_latency_ms=p95,
        violations=violations,
    )


