"""
Database layer for the multi-omics platform.

This package provides:
- SQLAlchemy models mapped to domain models
- Database connection and session management
- Alembic migration support
"""

from amprenta_rag.database.base import Base, get_db, get_engine, get_session_local, init_db

# Import models to register them with Base
from amprenta_rag.database import models  # noqa: F401

# Don't initialize SessionLocal on import - let it be lazy
# SessionLocal will be created when get_session_local() is first called

__all__ = [
    "Base",
    "SessionLocal",
    "get_db",
    "get_engine",
    "get_session_local",
    "init_db",
    "models",
]

