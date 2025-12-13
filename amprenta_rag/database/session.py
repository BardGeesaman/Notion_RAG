"""Database session management context manager."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Generator

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db


@contextmanager
def db_session() -> Generator[Session, None, None]:
    """
    Context manager for database sessions with automatic cleanup.

    Ensures proper rollback on errors and session closure to prevent connection leaks.
    """
    db_gen = None
    db = None

    try:
        db_gen = get_db()
        db = next(db_gen)

        # Ensure session is in a clean state before use
        try:
            db.rollback()
        except Exception:
            try:
                if db:
                    db.close()
            except Exception:
                pass
            db_gen = get_db()
            db = next(db_gen)

        yield db
        db.commit()

    except Exception as e:
        if db:
            try:
                db.rollback()
            except Exception:
                pass
        raise e

    finally:
        if db:
            try:
                db.close()
            except Exception:
                pass
        if db_gen:
            try:
                next(db_gen, None)
            except StopIteration:
                pass
            except Exception:
                pass

