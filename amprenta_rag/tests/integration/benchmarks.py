"""Performance benchmark utilities for integration tests."""

import json
import time
import logging
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, asdict, field
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResult:
    """Single benchmark measurement."""
    test_name: str
    endpoint: str
    method: str
    elapsed_ms: float
    threshold_ms: float
    passed: bool
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())


class BenchmarkTracker:
    """Track and report performance benchmarks."""
    
    # Thresholds per operation type (milliseconds)
    THRESHOLDS = {
        "GET": 150,       # Simple reads
        "POST": 200,      # Creates
        "PUT": 200,       # Updates
        "PATCH": 200,     # Partial updates
        "DELETE": 100,    # Deletes
        "batch": 500,     # Batch operations
        "compute": 3000,  # ML/compute operations
    }
    
    def __init__(self):
        self.results: list[BenchmarkResult] = []
    
    def _get_threshold(self, method: str, endpoint: str) -> float:
        """Get threshold for operation type."""
        # Check for special endpoints
        if "batch" in endpoint.lower():
            return self.THRESHOLDS["batch"]
        if any(x in endpoint.lower() for x in ["predict", "compute", "admet", "qsar"]):
            return self.THRESHOLDS["compute"]
        return self.THRESHOLDS.get(method.upper(), 200)
    
    def record(self, test_name: str, endpoint: str, method: str, elapsed_ms: float) -> BenchmarkResult:
        """Record a benchmark measurement."""
        threshold = self._get_threshold(method, endpoint)
        result = BenchmarkResult(
            test_name=test_name,
            endpoint=endpoint,
            method=method,
            elapsed_ms=round(elapsed_ms, 2),
            threshold_ms=threshold,
            passed=elapsed_ms <= threshold,
        )
        self.results.append(result)
        
        if not result.passed:
            logger.warning(
                f"SLOW: {test_name} - {method} {endpoint} took {elapsed_ms:.1f}ms "
                f"(threshold: {threshold}ms)"
            )
        
        return result
    
    def report(self) -> dict:
        """Generate benchmark report."""
        passed = [r for r in self.results if r.passed]
        failed = [r for r in self.results if not r.passed]
        
        return {
            "summary": {
                "total": len(self.results),
                "passed": len(passed),
                "failed": len(failed),
                "pass_rate": round(len(passed) / len(self.results) * 100, 1) if self.results else 0,
            },
            "slowest": [
                asdict(r) for r in sorted(self.results, key=lambda r: r.elapsed_ms, reverse=True)[:10]
            ],
            "failures": [asdict(r) for r in failed],
            "all_results": [asdict(r) for r in self.results],
        }
    
    def save(self, path: Path):
        """Save report to JSON file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.report(), f, indent=2)
        logger.info(f"Benchmark report saved to {path}")


# Global tracker instance for pytest
_tracker: Optional[BenchmarkTracker] = None


def get_tracker() -> BenchmarkTracker:
    """Get or create global benchmark tracker."""
    global _tracker
    if _tracker is None:
        _tracker = BenchmarkTracker()
    return _tracker


def reset_tracker():
    """Reset global tracker (for test isolation)."""
    global _tracker
    _tracker = BenchmarkTracker()
