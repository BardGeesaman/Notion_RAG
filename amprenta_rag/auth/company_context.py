"""Tenant company context integration for Postgres RLS."""

from __future__ import annotations

from sqlalchemy import text
from sqlalchemy.orm import Session

from amprenta_rag.models.auth import User


def set_company_context(request, db: Session, current_user: User) -> User:  # noqa: ARG001
    """Set the current company id into the DB session for RLS policies.

    Executes:
      SET LOCAL app.current_company_id = '{user.company_id}'

    Returns:
      current_user (for downstream use)
    """
    if not getattr(current_user, "company_id", None):
        return current_user

    # SET LOCAL is scoped to the current transaction; SQLAlchemy session keeps it within the request.
    db.execute(text(f"SET LOCAL app.current_company_id = '{current_user.company_id}'"))
    return current_user


__all__ = ["set_company_context"]


