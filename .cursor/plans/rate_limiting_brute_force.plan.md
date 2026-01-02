# Rate Limiting & Brute Force Protection

## Overview

Add comprehensive rate limiting and brute force protection to the FastAPI API to prevent abuse, credential stuffing, and denial-of-service attacks.

## Current State

- **No rate limiting middleware** - `slowapi` not in requirements.txt
- **Manual rate limiting** exists only for external APIs (KEGG, VEP)
- **No account lockout** - Failed password attempts not tracked
- **No connection timeout protection** - Slowloris attacks possible
- **Password-protected endpoints** - `/signatures/sign` requires password verification
- **Public endpoints** - `/share-links/{token}/validate` has no auth

## Sensitive Endpoints Requiring Protection

| Endpoint | Risk | Recommended Limit |
|----------|------|-------------------|
| `POST /signatures/sign` | Password brute force | 5/minute per IP |
| `GET /share-links/{token}/validate` | Token enumeration | 30/minute per IP |
| `POST /api/v1/admin/*` | Admin abuse | 10/minute per user |
| `POST /api/v1/compounds` | Resource exhaustion | 100/minute per user |
| `POST /api/v1/datasets/*/upload` | Upload abuse | 10/minute per user |
| `POST /extraction/ocr` | Compute-intensive | 5/minute per user |
| `POST /extraction/scrape` | SSRF amplification | 10/minute per user |

## Implementation Plan

### Batch 1: Core Rate Limiting Infrastructure (Implementor)

1. **Add slowapi to requirements.txt**
   ```
   slowapi>=0.1.9
   ```

2. **Create `amprenta_rag/api/rate_limit.py`**
   ```python
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
       return JSONResponse(
           status_code=429,
           content={
               "error": "rate_limit_exceeded",
               "detail": f"Rate limit exceeded: {exc.detail}",
               "retry_after": exc.retry_after,
           },
           headers={"Retry-After": str(exc.retry_after)},
       )

   # Key functions for different contexts
   def get_user_or_ip(request: Request) -> str:
       """Get authenticated user ID if available, otherwise IP.
       
       P1 FIX: Do NOT trust X-User-Id header directly (can be spoofed).
       Instead, check request.state.user set by auth dependency.
       """
       # Check if user was authenticated via get_current_user dependency
       user = getattr(request.state, "user", None)
       if user and hasattr(user, "id"):
           return f"user:{user.id}"
       return f"ip:{get_remote_address(request)}"
   
   
   def set_user_state(request: Request, user) -> None:
       """Set authenticated user on request state (call from get_current_user)."""
       request.state.user = user
   ```

3. **Integrate into `amprenta_rag/api/main.py`**
   ```python
   from slowapi import _rate_limit_exceeded_handler
   from slowapi.errors import RateLimitExceeded
   from amprenta_rag.api.rate_limit import limiter, rate_limit_exceeded_handler

   app.state.limiter = limiter
   app.add_exception_handler(RateLimitExceeded, rate_limit_exceeded_handler)
   ```

### Batch 2: Protect Sensitive Endpoints (Implementor)

1. **Signature endpoint** (`routers/signatures.py`)
   ```python
   from amprenta_rag.api.rate_limit import limiter

   @router.post("/sign", response_model=SignResponse, status_code=201)
   @limiter.limit("5/minute")
   async def sign_document(request: Request, ...):
   ```

2. **Share link validation** (`routers/share_links.py`)
   ```python
   @router.get("/{token}/validate", response_model=ShareLinkValidateResponse)
   @limiter.limit("30/minute")
   async def validate_token(request: Request, ...):
   ```

3. **Admin endpoints** (`routers/admin.py`)
   ```python
   @router.post("/...")
   @limiter.limit("10/minute", key_func=get_user_or_ip)
   async def admin_operation(request: Request, ...):
   ```

4. **Upload/compute-intensive endpoints**
   - OCR: 5/minute
   - Scrape: 10/minute
   - Dataset upload: 10/minute

### Batch 3: Account Lockout (Implementor)

**P1 FIX: Use database table for lockout persistence (survives restarts)**

1. **Create Alembic migration for lockout table**
   ```python
   # alembic/versions/xxx_add_auth_lockout_table.py
   def upgrade():
       op.create_table(
           'auth_lockout_attempts',
           sa.Column('id', sa.Integer(), primary_key=True),
           sa.Column('identifier', sa.String(255), nullable=False, index=True),
           sa.Column('attempt_time', sa.DateTime(timezone=True), nullable=False),
           sa.Column('created_at', sa.DateTime(timezone=True), server_default=sa.func.now()),
       )
       op.create_index('ix_lockout_identifier_time', 'auth_lockout_attempts', ['identifier', 'attempt_time'])
   ```

2. **Create `amprenta_rag/auth/lockout.py`**
   ```python
   """Account lockout tracking for brute force protection.
   
   Uses database for persistence (P1 fix: survives restarts).
   """
   from datetime import datetime, timedelta, timezone
   from typing import Optional, Tuple
   import logging
   
   from sqlalchemy import Column, Integer, String, DateTime, func, and_
   from sqlalchemy.orm import Session
   
   from amprenta_rag.database.base import Base
   from amprenta_rag.database.session import db_session
   
   logger = logging.getLogger(__name__)
   
   MAX_FAILED_ATTEMPTS = 5
   LOCKOUT_DURATION = timedelta(minutes=15)
   ATTEMPT_WINDOW = timedelta(minutes=5)
   
   
   class AuthLockoutAttempt(Base):
       """Track failed authentication attempts."""
       __tablename__ = 'auth_lockout_attempts'
       
       id = Column(Integer, primary_key=True)
       identifier = Column(String(255), nullable=False, index=True)
       attempt_time = Column(DateTime(timezone=True), nullable=False)
   
   
   def record_failed_attempt(identifier: str, db: Optional[Session] = None) -> None:
       """Record a failed authentication attempt."""
       now = datetime.now(timezone.utc)
       
       def _record(session: Session):
           attempt = AuthLockoutAttempt(identifier=identifier, attempt_time=now)
           session.add(attempt)
           session.commit()
           
           # Clean old attempts
           cutoff = now - ATTEMPT_WINDOW
           session.query(AuthLockoutAttempt).filter(
               AuthLockoutAttempt.identifier == identifier,
               AuthLockoutAttempt.attempt_time < cutoff
           ).delete()
           session.commit()
           
           # Log count
           count = session.query(AuthLockoutAttempt).filter(
               AuthLockoutAttempt.identifier == identifier,
               AuthLockoutAttempt.attempt_time >= cutoff
           ).count()
           logger.warning(f"Failed auth attempt for {identifier}: {count}/{MAX_FAILED_ATTEMPTS}")
       
       if db:
           _record(db)
       else:
           with db_session() as session:
               _record(session)
   
   
   def is_locked_out(identifier: str, db: Optional[Session] = None) -> Tuple[bool, Optional[int]]:
       """Check if identifier is locked out.
       
       Returns:
           (is_locked, seconds_remaining)
       """
       now = datetime.now(timezone.utc)
       cutoff = now - ATTEMPT_WINDOW
       
       def _check(session: Session) -> Tuple[bool, Optional[int]]:
           count = session.query(AuthLockoutAttempt).filter(
               AuthLockoutAttempt.identifier == identifier,
               AuthLockoutAttempt.attempt_time >= cutoff
           ).count()
           
           if count >= MAX_FAILED_ATTEMPTS:
               # Get last attempt time
               last = session.query(func.max(AuthLockoutAttempt.attempt_time)).filter(
                   AuthLockoutAttempt.identifier == identifier
               ).scalar()
               
               if last:
                   lockout_end = last + LOCKOUT_DURATION
                   if now < lockout_end:
                       remaining = int((lockout_end - now).total_seconds())
                       return True, remaining
           
           return False, None
       
       if db:
           return _check(db)
       else:
           with db_session() as session:
               return _check(session)
   
   
   def clear_failed_attempts(identifier: str, db: Optional[Session] = None) -> None:
       """Clear failed attempts after successful auth."""
       def _clear(session: Session):
           session.query(AuthLockoutAttempt).filter(
               AuthLockoutAttempt.identifier == identifier
           ).delete()
           session.commit()
       
       if db:
           _clear(db)
       else:
           with db_session() as session:
               _clear(session)
   ```

2. **Integrate into signature endpoint**
   ```python
   from amprenta_rag.auth.lockout import record_failed_attempt, is_locked_out, clear_failed_attempts

   @router.post("/sign", ...)
   async def sign_document(...):
       identifier = f"sign:{request.client.host}"
       
       locked, remaining = is_locked_out(identifier)
       if locked:
           raise HTTPException(
               status_code=429,
               detail=f"Account locked due to too many failed attempts. Try again in {remaining} seconds.",
               headers={"Retry-After": str(remaining)},
           )
       
       signature = create_signature(...)
       
       if not signature:
           record_failed_attempt(identifier)
           raise HTTPException(status_code=401, detail="Authentication failed")
       
       clear_failed_attempts(identifier)
       return SignResponse(...)
   ```

### Batch 4: Slowloris Protection (Implementor)

1. **Add request timeout middleware** (`api/middleware.py`)
   ```python
   """Custom middleware for security."""
   import asyncio
   import os
   from fastapi import Request
   from starlette.middleware.base import BaseHTTPMiddleware
   from starlette.responses import JSONResponse
   
   REQUEST_TIMEOUT = int(os.getenv("REQUEST_TIMEOUT_SECONDS", "60"))
   
   class TimeoutMiddleware(BaseHTTPMiddleware):
       """Timeout long-running requests to prevent Slowloris attacks."""
       
       def __init__(self, app, timeout: int = REQUEST_TIMEOUT):
           super().__init__(app)
           self.timeout = timeout
       
       async def dispatch(self, request: Request, call_next):
           try:
               return await asyncio.wait_for(
                   call_next(request),
                   timeout=self.timeout
               )
           except asyncio.TimeoutError:
               return JSONResponse(
                   status_code=504,
                   content={"error": "request_timeout", "detail": "Request timed out"}
               )
   ```

2. **Configure in main.py**
   ```python
   from amprenta_rag.api.middleware import TimeoutMiddleware
   
   # Default 60s, configurable via REQUEST_TIMEOUT_SECONDS env var
   app.add_middleware(TimeoutMiddleware)
   ```

### Batch 5: Tests + Documentation (Implementor)

1. **Create `amprenta_rag/tests/api/test_rate_limiting.py`**
   - Test rate limit exceeded returns 429
   - Test rate limit headers present
   - Test lockout after failed attempts
   - Test lockout expiry
   - Test successful auth clears lockout

2. **Create `amprenta_rag/tests/auth/test_lockout.py`**
   - Unit tests for lockout module

3. **Update `docs/ROADMAP.md`**
   - Mark "Rate Limiting & Brute Force Protection" as COMPLETE

4. **Update `.env.example`**
   - Add `REDIS_URL` for production rate limiting

## Key Files

| File | Action |
|------|--------|
| `requirements.txt` | ADD slowapi |
| `amprenta_rag/api/rate_limit.py` | CREATE |
| `amprenta_rag/api/middleware.py` | CREATE |
| `amprenta_rag/api/main.py` | UPDATE |
| `amprenta_rag/auth/lockout.py` | CREATE |
| `amprenta_rag/api/routers/signatures.py` | UPDATE |
| `amprenta_rag/api/routers/share_links.py` | UPDATE |
| `amprenta_rag/api/routers/extraction.py` | UPDATE |
| `amprenta_rag/tests/api/test_rate_limiting.py` | CREATE |
| `amprenta_rag/tests/auth/test_lockout.py` | CREATE |
| `docs/ROADMAP.md` | UPDATE |
| `.env.example` | UPDATE |

## Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `REDIS_URL` | None | Redis URL for distributed rate limiting (memory if not set) |
| `RATE_LIMIT_DEFAULT` | 200/minute | Default rate limit for all endpoints |
| `LOCKOUT_MAX_ATTEMPTS` | 5 | Failed attempts before lockout |
| `LOCKOUT_DURATION_MINUTES` | 15 | Lockout duration |
| `REQUEST_TIMEOUT_SECONDS` | 60 | Maximum request duration |

## Effort Estimate

- **Batch 1**: 2 hours (core infrastructure)
- **Batch 2**: 2 hours (endpoint protection)
- **Batch 3**: 2 hours (lockout system)
- **Batch 4**: 1 hour (timeout middleware)
- **Batch 5**: 2 hours (tests + docs)
- **Total**: ~9 hours (1-2 days)

## P2 Observations (Future)

- Redis backend for distributed rate limiting
- Per-endpoint configurable limits via config file
- Rate limit dashboard in admin UI
- Webhook alerts for repeated lockouts
- IP reputation scoring integration
- IP-based lockout may affect legitimate users behind corporate NAT (acceptable for MVP)
- Metrics/alerting for repeated lockout attempts

