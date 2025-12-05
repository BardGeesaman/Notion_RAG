"""
Production utilities for the Amprenta RAG system.

This package provides:
- Error handling and recovery (retry, circuit breakers, graceful degradation)
- Configuration validation
- Performance profiling and monitoring
- Health checks
"""

from amprenta_rag.utils.config_validation import (
    ConfigValidator,
    ValidationError,
    ValidationResult,
    validate_configuration,
)
from amprenta_rag.utils.error_handling import (
    CircuitBreaker,
    CircuitBreakerOpen,
    ErrorTracker,
    RetryConfig,
    RetryableError,
    get_error_tracker,
    graceful_degradation,
    retry_with_backoff,
)
from amprenta_rag.utils.health_check import (
    HealthCheckResult,
    HealthChecker,
    check_system_health,
)
from amprenta_rag.utils.performance import (
    PerformanceMetrics,
    PerformanceTimer,
    get_performance_metrics,
    timed,
    timer,
)

__all__ = [
    # Error handling
    "RetryConfig",
    "retry_with_backoff",
    "CircuitBreaker",
    "CircuitBreakerOpen",
    "RetryableError",
    "graceful_degradation",
    "ErrorTracker",
    "get_error_tracker",
    # Configuration validation
    "ConfigValidator",
    "ValidationError",
    "ValidationResult",
    "validate_configuration",
    # Health checks
    "HealthChecker",
    "HealthCheckResult",
    "check_system_health",
    # Performance
    "PerformanceTimer",
    "PerformanceMetrics",
    "get_performance_metrics",
    "timed",
    "timer",
]

