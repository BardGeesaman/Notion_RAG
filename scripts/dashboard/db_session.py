"""Database session management for the Streamlit dashboard."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Generator

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db


@contextmanager
def db_session() -> Generator[Session, None, None]:
    """
    Context manager for database sessions with automatic cleanup.

    Ensures proper rollback on errors and session closure.
    This is critical for preventing connection leaks and database locks.

    Features:
    - Automatic rollback on errors
    - Explicit commit on success
    - Guaranteed connection closure in finally block
    - Handles invalid transaction states

    Yields:
        Session: SQLAlchemy database session

    Usage:
        with db_session() as db:
            programs = db.query(Program).all()
    """
    db_gen = None
    db = None

    try:
        db_gen = get_db()
        db = next(db_gen)

        # Ensure session is in a clean state before use
        try:
            db.rollback()  # Clear any pending transaction
        except Exception:
            # If rollback fails, the session might be in a bad state
            # Close it and get a fresh one
            try:
                if db:
                    db.close()
            except Exception:
                pass
            # Get a new session
            db_gen = get_db()
            db = next(db_gen)

        # Yield the session for use
        yield db

        # Commit successful transactions
        db.commit()

    except Exception as e:
        # Rollback on any error to prevent leaving transaction in invalid state
        if db:
            try:
                db.rollback()
            except Exception:
                # If rollback fails, the connection is likely broken
                pass
        raise e

    finally:
        # CRITICAL: Always close the session to release the connection
        # This prevents connection leaks and database locks
        if db:
            try:
                db.close()
            except Exception:
                # Even if close fails, we tried
                pass

        # Also ensure the generator is properly exhausted
        # This triggers the finally block in get_db() which also closes the session
        if db_gen:
            try:
                next(db_gen, None)  # Exhaust generator to trigger its finally block
            except StopIteration:
                pass
            except Exception:
                # If generator is already exhausted or broken, that's okay
                pass
