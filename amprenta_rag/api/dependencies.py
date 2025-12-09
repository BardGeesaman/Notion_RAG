"""
FastAPI dependencies for database sessions and authentication.

Provides reusable dependencies for API routes.
"""

from typing import Generator

from sqlalchemy.orm import Session

from amprenta_rag.database.base import get_db


def get_database_session() -> Generator[Session, None, None]:
    """
    Dependency to get a database session.
    
    Yields:
        Database session
    """
    yield from get_db()

