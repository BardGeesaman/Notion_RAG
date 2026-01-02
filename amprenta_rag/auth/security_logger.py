"""Security event logging for audit and monitoring.

Provides structured logging for security-relevant events including:
- Authentication successes and failures
- Authorization failures
- Rate limit violations
- Account lockouts
- Input validation failures
- Suspicious activity patterns
"""

import logging
from datetime import datetime, timezone
from enum import Enum
from typing import Optional
from uuid import UUID

# Dedicated security logger
security_logger = logging.getLogger("security")
security_logger.setLevel(logging.INFO)


class SecurityEvent(str, Enum):
    """Security event types for classification."""
    
    # Authentication events
    AUTH_SUCCESS = "auth_success"
    AUTH_FAILURE = "auth_failure"
    AUTH_LOGOUT = "auth_logout"
    
    # Authorization events
    AUTHZ_FAILURE = "authz_failure"
    PERMISSION_DENIED = "permission_denied"
    
    # Rate limiting
    RATE_LIMIT_HIT = "rate_limit_hit"
    RATE_LIMIT_EXCEEDED = "rate_limit_exceeded"
    
    # Account security
    ACCOUNT_LOCKED = "account_locked"
    ACCOUNT_UNLOCKED = "account_unlocked"
    PASSWORD_CHANGED = "password_changed"
    
    # Input validation
    VALIDATION_FAILURE = "validation_failure"
    INJECTION_ATTEMPT = "injection_attempt"
    
    # Suspicious activity
    SUSPICIOUS_ACTIVITY = "suspicious_activity"
    ENUMERATION_ATTEMPT = "enumeration_attempt"
    
    # Data access
    SENSITIVE_DATA_ACCESS = "sensitive_data_access"
    BULK_DATA_EXPORT = "bulk_data_export"


def log_security_event(
    event: SecurityEvent,
    user_id: Optional[UUID] = None,
    ip_address: Optional[str] = None,
    endpoint: Optional[str] = None,
    details: Optional[dict] = None,
    severity: str = "INFO"
) -> None:
    """
    Log a security event with structured data.
    
    Args:
        event: Type of security event
        user_id: User ID if authenticated
        ip_address: Client IP address
        endpoint: API endpoint accessed
        details: Additional event details
        severity: Log level (INFO, WARNING, ERROR)
    """
    log_data = {
        "event_type": event.value,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "user_id": str(user_id) if user_id else None,
        "ip_address": ip_address,
        "endpoint": endpoint,
        "details": details or {},
    }
    
    # Format as structured log message
    message = (
        f"[SECURITY] {event.value} | "
        f"user={user_id or 'anonymous'} | "
        f"ip={ip_address or 'unknown'} | "
        f"endpoint={endpoint or 'N/A'}"
    )
    
    if details:
        message += f" | details={details}"
    
    # Log at appropriate level
    log_level = getattr(logging, severity.upper(), logging.INFO)
    security_logger.log(log_level, message, extra={"security_data": log_data})


def log_auth_failure(
    username: Optional[str] = None,
    ip_address: Optional[str] = None,
    reason: str = "invalid_credentials"
) -> None:
    """Log authentication failure."""
    log_security_event(
        SecurityEvent.AUTH_FAILURE,
        ip_address=ip_address,
        details={"username": username, "reason": reason},
        severity="WARNING"
    )


def log_auth_success(
    user_id: UUID,
    ip_address: Optional[str] = None
) -> None:
    """Log successful authentication."""
    log_security_event(
        SecurityEvent.AUTH_SUCCESS,
        user_id=user_id,
        ip_address=ip_address,
        severity="INFO"
    )


def log_authz_failure(
    user_id: Optional[UUID],
    resource: str,
    action: str,
    ip_address: Optional[str] = None
) -> None:
    """Log authorization failure (403)."""
    log_security_event(
        SecurityEvent.AUTHZ_FAILURE,
        user_id=user_id,
        ip_address=ip_address,
        details={"resource": resource, "action": action},
        severity="WARNING"
    )


def log_rate_limit(
    user_id: Optional[UUID],
    endpoint: str,
    ip_address: Optional[str] = None,
    limit: Optional[str] = None
) -> None:
    """Log rate limit violation."""
    log_security_event(
        SecurityEvent.RATE_LIMIT_HIT,
        user_id=user_id,
        ip_address=ip_address,
        endpoint=endpoint,
        details={"limit": limit},
        severity="WARNING"
    )


def log_account_locked(
    user_id: UUID,
    ip_address: Optional[str] = None,
    failed_attempts: int = 0
) -> None:
    """Log account lockout."""
    log_security_event(
        SecurityEvent.ACCOUNT_LOCKED,
        user_id=user_id,
        ip_address=ip_address,
        details={"failed_attempts": failed_attempts},
        severity="WARNING"
    )


def log_validation_failure(
    endpoint: str,
    field: str,
    value_type: str,
    ip_address: Optional[str] = None
) -> None:
    """Log input validation failure."""
    log_security_event(
        SecurityEvent.VALIDATION_FAILURE,
        ip_address=ip_address,
        endpoint=endpoint,
        details={"field": field, "value_type": value_type},
        severity="INFO"
    )


def log_suspicious_activity(
    description: str,
    user_id: Optional[UUID] = None,
    ip_address: Optional[str] = None,
    details: Optional[dict] = None
) -> None:
    """Log suspicious activity for review."""
    log_security_event(
        SecurityEvent.SUSPICIOUS_ACTIVITY,
        user_id=user_id,
        ip_address=ip_address,
        details={"description": description, **(details or {})},
        severity="WARNING"
    )
