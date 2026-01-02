"""Rate limiting configuration and utilities."""
from slowapi import Limiter
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
from fastapi import Request
from fastapi.responses import JSONResponse
import os

# Use Redis in production, memory in development
REDIS_URL = os.getenv("REDIS_URL")

limiter = Limiter(
    key_func=get_remote_address,
    storage_uri=REDIS_URL if REDIS_URL else "memory://",
    default_limits=["200/minute"],  # Global default
)


def rate_limit_exceeded_handler(request: Request, exc: RateLimitExceeded):
    """Custom handler for rate limit exceeded."""
    return JSONResponse(
        status_code=429,
        content={
            "error": "rate_limit_exceeded",
            "detail": f"Rate limit exceeded: {exc.detail}",
            "retry_after": exc.retry_after,
        },
        headers={"Retry-After": str(exc.retry_after)},
    )


def get_user_or_ip(request: Request) -> str:
    """Get authenticated user ID if available, otherwise IP.
    
    P1 FIX: Do NOT trust X-User-Id header directly (can be spoofed).
    Instead, check request.state.user set by auth dependency.
    """
    user = getattr(request.state, "user", None)
    if user and hasattr(user, "id"):
        return f"user:{user.id}"
    return f"ip:{get_remote_address(request)}"


def set_user_state(request: Request, user) -> None:
    """Set authenticated user on request state (call from get_current_user)."""
    request.state.user = user
