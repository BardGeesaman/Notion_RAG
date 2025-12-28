"""
Unit tests for RateLimiter.

Tests timing, thread safety, and NCBI API key detection.
"""

from __future__ import annotations

import os
import threading
import time
from unittest.mock import patch

import pytest

from amprenta_rag.ingestion.papers.rate_limiter import RateLimiter


class TestRateLimiter:
    """Tests for RateLimiter class."""

    def test_rate_limiter_initialization(self):
        """Test rate limiter initializes with correct parameters."""
        limiter = RateLimiter(requests_per_second=5.0)

        assert limiter.requests_per_second == 5.0
        assert limiter.min_interval == pytest.approx(0.2, abs=0.01)

    def test_rate_limiter_enforces_delay(self):
        """Test rate limiter enforces minimum delay between calls."""
        limiter = RateLimiter(requests_per_second=10.0)  # 0.1s interval

        @limiter
        def test_func():
            return time.time()

        # First call
        start = time.time()
        t1 = test_func()

        # Second call should be delayed
        t2 = test_func()
        elapsed = t2 - t1

        # Should be approximately 0.1 seconds
        assert elapsed >= 0.09  # Allow small margin
        assert elapsed < 0.15

    def test_rate_limiter_multiple_calls(self):
        """Test rate limiter works correctly across multiple calls."""
        limiter = RateLimiter(requests_per_second=20.0)  # 0.05s interval

        @limiter
        def test_func():
            pass

        start = time.time()

        # Make 5 calls
        for _ in range(5):
            test_func()

        elapsed = time.time() - start

        # Should take approximately 4 * 0.05 = 0.2 seconds (4 intervals)
        assert elapsed >= 0.18
        assert elapsed < 0.3

    def test_rate_limiter_decorator_preserves_return_value(self):
        """Test decorator preserves function return value."""
        limiter = RateLimiter(requests_per_second=100.0)

        @limiter
        def get_value():
            return 42

        assert get_value() == 42

    def test_rate_limiter_decorator_preserves_arguments(self):
        """Test decorator preserves function arguments."""
        limiter = RateLimiter(requests_per_second=100.0)

        @limiter
        def add(a, b):
            return a + b

        assert add(3, 4) == 7
        assert add(10, 20) == 30

    def test_rate_limiter_thread_safety(self):
        """Test rate limiter is thread-safe."""
        limiter = RateLimiter(requests_per_second=20.0)
        call_times = []
        lock = threading.Lock()

        @limiter
        def test_func():
            with lock:
                call_times.append(time.time())

        # Create multiple threads
        threads = []
        for _ in range(5):
            t = threading.Thread(target=test_func)
            threads.append(t)

        # Start all threads at once
        start = time.time()
        for t in threads:
            t.start()

        # Wait for all threads
        for t in threads:
            t.join()

        # Verify timing
        call_times.sort()
        for i in range(1, len(call_times)):
            interval = call_times[i] - call_times[i - 1]
            # Each call should be separated by minimum interval
            assert interval >= 0.04  # Allow small margin for 0.05s

    def test_for_ncbi_without_api_key(self, monkeypatch):
        """Test for_ncbi() creates 3 req/sec limiter without API key."""
        # Ensure no API key
        monkeypatch.delenv("NCBI_API_KEY", raising=False)

        limiter = RateLimiter.for_ncbi()

        assert limiter.requests_per_second == 3.0
        assert limiter.min_interval == pytest.approx(1.0 / 3.0, abs=0.01)

    def test_for_ncbi_with_api_key(self, monkeypatch):
        """Test for_ncbi() creates 10 req/sec limiter with API key."""
        monkeypatch.setenv("NCBI_API_KEY", "test_api_key_123")

        limiter = RateLimiter.for_ncbi()

        assert limiter.requests_per_second == 10.0
        assert limiter.min_interval == pytest.approx(0.1, abs=0.01)

    def test_rate_limiter_with_fast_function(self):
        """Test rate limiter with function that executes quickly."""
        limiter = RateLimiter(requests_per_second=10.0)
        call_count = 0

        @limiter
        def fast_func():
            nonlocal call_count
            call_count += 1
            return call_count

        start = time.time()

        # Make 3 calls
        results = [fast_func() for _ in range(3)]

        elapsed = time.time() - start

        # Should take approximately 2 * 0.1 = 0.2 seconds (2 intervals)
        assert elapsed >= 0.18
        assert elapsed < 0.3
        assert results == [1, 2, 3]

    def test_rate_limiter_with_slow_function(self):
        """Test rate limiter with function that takes time to execute."""
        limiter = RateLimiter(requests_per_second=10.0)

        @limiter
        def slow_func():
            time.sleep(0.05)  # Function takes 50ms
            return True

        start = time.time()

        # Make 2 calls
        slow_func()
        slow_func()

        elapsed = time.time() - start

        # First call: 50ms (function)
        # Second call: 50ms wait (to reach 100ms interval) + 50ms (function) = 150ms total
        assert elapsed >= 0.14
        assert elapsed < 0.20

