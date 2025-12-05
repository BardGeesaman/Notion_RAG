"""
Enhanced error handling and recovery utilities for production use.

Provides retry logic, circuit breakers, graceful degradation,
and structured error reporting.
"""

from __future__ import annotations

import functools
import time
from typing import Any, Callable, Optional, Type, TypeVar

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

T = TypeVar("T")


class RetryableError(Exception):
    """Base exception for errors that should trigger retries."""
    pass


class CircuitBreakerOpen(Exception):
    """Raised when circuit breaker is open."""
    pass


class RetryConfig:
    """Configuration for retry logic."""
    
    def __init__(
        self,
        max_attempts: int = 3,
        initial_delay: float = 1.0,
        max_delay: float = 60.0,
        exponential_base: float = 2.0,
        retryable_exceptions: tuple[Type[Exception], ...] = (Exception,),
    ):
        self.max_attempts = max_attempts
        self.initial_delay = initial_delay
        self.max_delay = max_delay
        self.exponential_base = exponential_base
        self.retryable_exceptions = retryable_exceptions


def retry_with_backoff(
    config: Optional[RetryConfig] = None,
    log_errors: bool = True,
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to retry a function with exponential backoff.
    
    Args:
        config: Retry configuration (uses defaults if None)
        log_errors: Whether to log each retry attempt
        
    Returns:
        Decorated function with retry logic
    """
    if config is None:
        config = RetryConfig()
    
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            delay = config.initial_delay
            last_exception: Optional[Exception] = None
            
            for attempt in range(1, config.max_attempts + 1):
                try:
                    return func(*args, **kwargs)
                except config.retryable_exceptions as e:
                    last_exception = e
                    
                    if attempt < config.max_attempts:
                        if log_errors:
                            logger.warning(
                                "[ERROR-HANDLER] Attempt %d/%d failed for %s: %r. "
                                "Retrying in %.1fs...",
                                attempt,
                                config.max_attempts,
                                func.__name__,
                                e,
                                delay,
                            )
                        time.sleep(delay)
                        delay = min(delay * config.exponential_base, config.max_delay)
                    else:
                        if log_errors:
                            logger.error(
                                "[ERROR-HANDLER] All %d attempts failed for %s: %r",
                                config.max_attempts,
                                func.__name__,
                                e,
                            )
            
            if last_exception:
                raise last_exception
            
            raise RuntimeError("Retry logic failed unexpectedly")
        
        return wrapper
    
    return decorator


class CircuitBreaker:
    """
    Circuit breaker pattern for external service calls.
    
    Opens circuit after failure threshold, prevents cascading failures.
    """
    
    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: float = 60.0,
        expected_exception: Type[Exception] = Exception,
    ):
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.expected_exception = expected_exception
        self.failure_count = 0
        self.last_failure_time: Optional[float] = None
        self.state = "closed"  # closed, open, half_open
    
    def call(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        """
        Execute function with circuit breaker protection.
        
        Args:
            func: Function to execute
            *args: Positional arguments
            **kwargs: Keyword arguments
            
        Returns:
            Function result
            
        Raises:
            CircuitBreakerOpen: If circuit is open
            Exception: If function call fails
        """
        # Check circuit state
        if self.state == "open":
            if self.last_failure_time:
                elapsed = time.time() - self.last_failure_time
                if elapsed >= self.recovery_timeout:
                    self.state = "half_open"
                    logger.info(
                        "[CIRCUIT-BREAKER] Moving %s to half-open state",
                        func.__name__,
                    )
                else:
                    raise CircuitBreakerOpen(
                        f"Circuit breaker is open for {func.__name__}. "
                        f"Wait {self.recovery_timeout - elapsed:.1f}s before retrying."
                    )
        
        # Attempt call
        try:
            result = func(*args, **kwargs)
            # Success: reset failure count
            if self.state == "half_open":
                logger.info(
                    "[CIRCUIT-BREAKER] %s recovered, closing circuit",
                    func.__name__,
                )
                self.state = "closed"
            self.failure_count = 0
            return result
        except self.expected_exception as e:
            self.failure_count += 1
            self.last_failure_time = time.time()
            
            if self.failure_count >= self.failure_threshold:
                self.state = "open"
                logger.error(
                    "[CIRCUIT-BREAKER] Circuit breaker opened for %s after %d failures",
                    func.__name__,
                    self.failure_count,
                )
            
            raise


def graceful_degradation(
    fallback_value: Any = None,
    fallback_func: Optional[Callable] = None,
    log_warning: bool = True,
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to provide graceful degradation on failure.
    
    Args:
        fallback_value: Value to return on failure
        fallback_func: Function to call on failure (takes same args)
        log_warning: Whether to log warnings
        
    Returns:
        Decorated function with graceful degradation
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if log_warning:
                    logger.warning(
                        "[GRACEFUL-DEGRADE] %s failed: %r. Using fallback.",
                        func.__name__,
                        e,
                    )
                
                if fallback_func:
                    try:
                        return fallback_func(*args, **kwargs)
                    except Exception as fallback_error:
                        logger.error(
                            "[GRACEFUL-DEGRADE] Fallback for %s also failed: %r",
                            func.__name__,
                            fallback_error,
                        )
                
                return fallback_value
        
        return wrapper
    
    return decorator


class ErrorTracker:
    """Track and report errors for monitoring."""
    
    def __init__(self, max_errors: int = 100):
        self.max_errors = max_errors
        self.errors: list[dict[str, Any]] = []
    
    def record_error(
        self,
        operation: str,
        error: Exception,
        context: Optional[dict[str, Any]] = None,
    ):
        """Record an error for tracking."""
        error_entry = {
            "timestamp": time.time(),
            "operation": operation,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "context": context or {},
        }
        
        self.errors.append(error_entry)
        
        # Keep only recent errors
        if len(self.errors) > self.max_errors:
            self.errors.pop(0)
    
    def get_error_summary(self) -> dict[str, Any]:
        """Get summary of recent errors."""
        if not self.errors:
            return {"total_errors": 0}
        
        error_types = {}
        for error in self.errors:
            error_type = error["error_type"]
            error_types[error_type] = error_types.get(error_type, 0) + 1
        
        return {
            "total_errors": len(self.errors),
            "error_types": error_types,
            "recent_errors": self.errors[-10:],  # Last 10 errors
        }


# Global error tracker instance
_global_error_tracker = ErrorTracker()


def get_error_tracker() -> ErrorTracker:
    """Get the global error tracker instance."""
    return _global_error_tracker

