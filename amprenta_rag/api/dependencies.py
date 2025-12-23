"""
FastAPI dependencies for database sessions and authentication.

Provides reusable dependencies for API routes.
"""

from typing import Generator, Optional
from uuid import UUID

from fastapi import Depends, HTTPException, Request
from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db
from amprenta_rag.auth.company_context import set_company_context
from amprenta_rag.database.models import User


def get_database_session() -> Generator[Session, None, None]:
    """
    Dependency to get a database session.

    Yields:
        Database session
    """
    yield from get_db()


def get_optional_user_id() -> Optional[UUID]:
    """Return user_id from auth context if available, else None.

    MVP implementation: returns None (allows anonymous subscriptions).
    When auth is implemented, extract user_id from request headers/session.
    """
    return None


def get_current_user(request: Request, db: Session = Depends(get_database_session)) -> User:
    """MVP current user resolution from request headers.

    Expected header:
      - X-User-Id: <uuid>

    This is a minimal implementation to support opt-in multi-tenancy context.
    """
    raw = request.headers.get("x-user-id") or request.headers.get("X-User-Id")
    if not raw:
        raise HTTPException(status_code=401, detail="Missing X-User-Id header")
    try:
        uid = UUID(str(raw))
    except Exception:
        raise HTTPException(status_code=401, detail="Invalid X-User-Id header")
    user = db.query(User).filter(User.id == uid).first()
    if not user:
        raise HTTPException(status_code=401, detail="User not found")
    return user


def get_current_user_with_company(
    request: Request,
    db: Session = Depends(get_database_session),
    current_user: User = Depends(get_current_user),
) -> User:
    """Return current user and set DB company context for RLS policies."""
    return set_company_context(request, db, current_user)

