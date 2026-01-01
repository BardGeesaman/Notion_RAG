"""Query timing utilities for performance monitoring."""

import time
import logging
from contextlib import contextmanager
from typing import Optional

logger = logging.getLogger(__name__)

SLOW_QUERY_THRESHOLD_MS = 100  # Log queries slower than 100ms


@contextmanager
def timed_query(operation: str, threshold_ms: Optional[float] = None):
    """Context manager to time and log slow database operations.
    
    Args:
        operation: Description of the operation being timed
        threshold_ms: Custom threshold in milliseconds (default: 100ms)
    
    Usage:
        with timed_query("fetch_programs"):
            programs = db.query(Program).all()
    """
    threshold = threshold_ms or SLOW_QUERY_THRESHOLD_MS
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed_ms = (time.perf_counter() - start) * 1000
        if elapsed_ms > threshold:
            logger.warning(
                f"SLOW QUERY: {operation} took {elapsed_ms:.1f}ms (threshold: {threshold}ms)"
            )
        else:
            logger.debug(f"Query: {operation} took {elapsed_ms:.1f}ms")


@contextmanager
def timed_api_call(operation: str, threshold_ms: Optional[float] = None):
    """Context manager to time and log slow API calls.
    
    Args:
        operation: Description of the API call being timed
        threshold_ms: Custom threshold in milliseconds (default: 200ms for API calls)
    
    Usage:
        with timed_api_call("POST /api/v1/programs"):
            response = client.post("/api/v1/programs", json=data)
    """
    # API calls have higher default threshold since they include network latency
    threshold = threshold_ms or 200
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed_ms = (time.perf_counter() - start) * 1000
        if elapsed_ms > threshold:
            logger.warning(
                f"SLOW API CALL: {operation} took {elapsed_ms:.1f}ms (threshold: {threshold}ms)"
            )
        else:
            logger.debug(f"API Call: {operation} took {elapsed_ms:.1f}ms")


def log_query_metrics(operation: str, elapsed_ms: float, 
                     record_count: Optional[int] = None,
                     threshold_ms: Optional[float] = None):
    """Log query metrics with optional record count.
    
    Args:
        operation: Description of the operation
        elapsed_ms: Time taken in milliseconds
        record_count: Number of records returned (optional)
        threshold_ms: Custom threshold for slow query detection
    """
    threshold = threshold_ms or SLOW_QUERY_THRESHOLD_MS
    
    if record_count is not None:
        per_record_ms = elapsed_ms / record_count if record_count > 0 else elapsed_ms
        metrics_msg = f"{operation} took {elapsed_ms:.1f}ms for {record_count} records ({per_record_ms:.2f}ms/record)"
    else:
        metrics_msg = f"{operation} took {elapsed_ms:.1f}ms"
    
    if elapsed_ms > threshold:
        logger.warning(f"SLOW QUERY: {metrics_msg} (threshold: {threshold}ms)")
    else:
        logger.debug(f"Query: {metrics_msg}")


class QueryTimer:
    """Class-based query timer for more complex timing scenarios."""
    
    def __init__(self, operation: str, threshold_ms: Optional[float] = None):
        self.operation = operation
        self.threshold = threshold_ms or SLOW_QUERY_THRESHOLD_MS
        self.start_time = None
        self.end_time = None
    
    def start(self):
        """Start timing."""
        self.start_time = time.perf_counter()
        return self
    
    def stop(self):
        """Stop timing and log result."""
        if self.start_time is None:
            logger.error(f"QueryTimer for {self.operation} was never started")
            return
        
        self.end_time = time.perf_counter()
        elapsed_ms = (self.end_time - self.start_time) * 1000
        
        if elapsed_ms > self.threshold:
            logger.warning(
                f"SLOW QUERY: {self.operation} took {elapsed_ms:.1f}ms (threshold: {self.threshold}ms)"
            )
        else:
            logger.debug(f"Query: {self.operation} took {elapsed_ms:.1f}ms")
        
        return elapsed_ms
    
    def __enter__(self):
        return self.start()
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()


# Convenience functions for common thresholds
def timed_complex_query(operation: str):
    """Time complex queries with higher threshold (500ms)."""
    return timed_query(operation, threshold_ms=500)


def timed_simple_query(operation: str):
    """Time simple queries with lower threshold (50ms)."""
    return timed_query(operation, threshold_ms=50)


def timed_bulk_operation(operation: str):
    """Time bulk operations with much higher threshold (2000ms)."""
    return timed_query(operation, threshold_ms=2000)
