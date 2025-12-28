"""
Rate limiter for API calls to prevent exceeding rate limits.

Provides a decorator-based, thread-safe rate limiter that can be applied
to any function to enforce a maximum requests per second limit.
"""

from __future__ import annotations

import os
import time
from functools import wraps
from threading import Lock
from typing import Any, Callable, TypeVar

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Type variable for generic function signatures
F = TypeVar('F', bound=Callable[..., Any])


class RateLimiter:
    """
    Thread-safe rate limiter using token bucket algorithm.

    Enforces a maximum number of requests per second across all threads.
    Designed for API rate limiting (e.g., NCBI E-utilities: 3 req/sec without
    API key, 10 req/sec with key).

    Example:
        >>> # Create rate limiter for NCBI (auto-detects API key)
        >>> ncbi_limiter = RateLimiter.for_ncbi()
        >>>
        >>> @ncbi_limiter
        ... def fetch_paper(pmid):
        ...     return requests.get(f"https://eutils.ncbi.nlm.nih.gov/...")
        >>>
        >>> # Or create custom rate limiter
        >>> limiter = RateLimiter(requests_per_second=5)
        >>> @limiter
        ... def my_api_call():
        ...     pass
    """

    def __init__(self, requests_per_second: float = 3.0):
        """
        Initialize rate limiter.

        Args:
            requests_per_second: Maximum requests allowed per second
        """
        self.requests_per_second = requests_per_second
        self.min_interval = 1.0 / requests_per_second
        self._last_call_time = 0.0
        self._lock = Lock()
        
        logger.debug(
            "[RATE_LIMITER] Initialized with %.2f req/sec (%.3f sec interval)",
            requests_per_second,
            self.min_interval,
        )

    @classmethod
    def for_ncbi(cls) -> RateLimiter:
        """
        Create rate limiter configured for NCBI E-utilities.

        Checks for NCBI_API_KEY environment variable and sets appropriate rate:
        - With API key: 10 requests/second
        - Without API key: 3 requests/second

        Returns:
            RateLimiter instance configured for NCBI
        """
        api_key = os.getenv("NCBI_API_KEY")
        if api_key:
            rate = 10.0
            logger.info("[RATE_LIMITER] NCBI_API_KEY detected - using 10 req/sec")
        else:
            rate = 3.0
            logger.info("[RATE_LIMITER] No NCBI_API_KEY - using 3 req/sec")
        
        return cls(requests_per_second=rate)

    def __call__(self, func: F) -> F:
        """
        Decorator to apply rate limiting to a function.

        Args:
            func: Function to wrap with rate limiting

        Returns:
            Wrapped function that respects rate limits
        """
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            with self._lock:
                # Calculate time since last call
                current_time = time.time()
                time_since_last_call = current_time - self._last_call_time
                
                # If not enough time has passed, sleep for the remainder
                if time_since_last_call < self.min_interval:
                    sleep_time = self.min_interval - time_since_last_call
                    logger.debug(
                        "[RATE_LIMITER] Sleeping %.3f sec (%.2f req/sec limit)",
                        sleep_time,
                        self.requests_per_second,
                    )
                    time.sleep(sleep_time)
                
                # Update last call time
                self._last_call_time = time.time()
            
            # Execute the function outside the lock to allow other threads
            # to queue up while this function executes
            return func(*args, **kwargs)
        
        return wrapper  # type: ignore[return-value]


# Convenience instance for NCBI API calls
ncbi_rate_limiter = RateLimiter.for_ncbi()

