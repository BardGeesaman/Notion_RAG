"""
Performance profiling and monitoring utilities.

Provides timing, memory tracking, and performance metrics
for production monitoring.
"""

from __future__ import annotations

import functools
import time
from contextlib import contextmanager
from typing import Any, Callable, Dict, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class PerformanceTimer:
    """Track execution time for operations."""

    def __init__(self, operation_name: str, log_threshold: float = 1.0):
        self.operation_name = operation_name
        self.log_threshold = log_threshold
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Context-manager protocol args are unused; keep signature for compatibility.
        del exc_type, exc_val, exc_tb
        self.end_time = time.time()
        elapsed = self.elapsed_time()

        if elapsed >= self.log_threshold:
            logger.info(
                "[PERF] %s took %.2f seconds",
                self.operation_name,
                elapsed,
            )

        return False

    def elapsed_time(self) -> float:
        """Get elapsed time in seconds."""
        if self.start_time is None:
            return 0.0
        end = self.end_time if self.end_time is not None else time.time()
        return end - self.start_time


def timed(
    operation_name: Optional[str] = None,
    log_threshold: float = 1.0,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """
    Decorator to time function execution.

    Args:
        operation_name: Name for the operation (defaults to function name)
        log_threshold: Only log if execution time exceeds this (seconds)

    Returns:
        Decorated function with timing
    """
    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        name = operation_name or func.__name__

        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            with PerformanceTimer(name, log_threshold):
                return func(*args, **kwargs)

        return wrapper

    return decorator


@contextmanager
def timer(operation_name: str, log_threshold: float = 1.0):
    """
    Context manager for timing operations.

    Args:
        operation_name: Name of the operation
        log_threshold: Only log if execution time exceeds this (seconds)

    Example:
        with timer("database_query"):
            result = query_database()
    """
    with PerformanceTimer(operation_name, log_threshold):
        yield


class PerformanceMetrics:
    """Track and aggregate performance metrics."""

    def __init__(self):
        self.metrics: Dict[str, list[float]] = {}

    def record(self, operation: str, duration: float):
        """Record a performance metric."""
        if operation not in self.metrics:
            self.metrics[operation] = []
        self.metrics[operation].append(duration)

    def get_statistics(self, operation: str) -> Optional[Dict[str, float]]:
        """Get statistics for an operation."""
        if operation not in self.metrics or not self.metrics[operation]:
            return None

        durations = self.metrics[operation]
        return {
            "count": len(durations),
            "mean": sum(durations) / len(durations),
            "min": min(durations),
            "max": max(durations),
            "total": sum(durations),
        }

    def get_summary(self) -> Dict[str, Dict[str, float]]:
        """Get summary of all metrics."""
        summary = {}
        for operation in self.metrics:
            stats = self.get_statistics(operation)
            if stats:
                summary[operation] = stats
        return summary

    def reset(self):
        """Clear all metrics."""
        self.metrics.clear()


# Global performance metrics instance
_global_metrics = PerformanceMetrics()


def get_performance_metrics() -> PerformanceMetrics:
    """Get the global performance metrics instance."""
    return _global_metrics

