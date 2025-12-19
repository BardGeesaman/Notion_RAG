"""
Database base configuration and session management.

Provides SQLAlchemy base, engine, and session management for Postgres.
"""

from __future__ import annotations

from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool

from amprenta_rag.config import get_config

# SQLAlchemy declarative base for models
Base = declarative_base()

# Global engine and session factory (initialized on first use)
_engine = None
_SessionLocal = None


def get_engine():
    """
    Get or create the database engine.

    Returns:
        SQLAlchemy engine instance
    """
    global _engine

    if _engine is None:
        cfg = get_config()
        postgres_cfg = cfg.postgres

        # Build connection string from config
        # Use full URL if provided, otherwise construct from components
        if postgres_cfg.url:
            db_url = postgres_cfg.url
        else:
            db_url = (
                f"postgresql://{postgres_cfg.user}:{postgres_cfg.password}"
                f"@{postgres_cfg.host}:{postgres_cfg.port}/{postgres_cfg.db}"
            )

        _engine = create_engine(
            db_url,
            poolclass=NullPool,  # Use NullPool for now, can switch to connection pooling later
            echo=postgres_cfg.echo,  # Set to True for SQL logging
        )

    return _engine


def get_session_local():
    """
    Get or create the session factory.

    Returns:
        Session factory class
    """
    global _SessionLocal

    if _SessionLocal is None:
        _SessionLocal = sessionmaker(
            autocommit=False,
            autoflush=False,
            bind=get_engine(),
        )

    return _SessionLocal


def get_db() -> Generator:
    """
    Get a database session.

    Yields:
        Database session
    """
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db() -> None:
    """
    Initialize the database (create all tables).

    This should be called after all models are imported.
    """
    # Import all models here to ensure they're registered with Base
    from amprenta_rag.database import models  # noqa: F401

    Base.metadata.create_all(bind=get_engine())

