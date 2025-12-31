"""Async database base configuration and session management."""

from __future__ import annotations

import os
from typing import AsyncGenerator

from sqlalchemy.ext.asyncio import (
    AsyncEngine,
    AsyncSession,
    async_sessionmaker,
    create_async_engine,
)

from amprenta_rag.config import get_config

_async_engine: AsyncEngine | None = None
_AsyncSessionLocal: async_sessionmaker[AsyncSession] | None = None


def get_async_engine() -> AsyncEngine:
    """Get or create the async database engine."""
    global _async_engine

    if _async_engine is None:
        cfg = get_config()
        postgres_cfg = cfg.postgres

        # Build async connection string
        if postgres_cfg.url:
            db_url = postgres_cfg.url.replace("postgresql://", "postgresql+asyncpg://")
        else:
            db_url = (
                f"postgresql+asyncpg://{postgres_cfg.user}:{postgres_cfg.password}"
                f"@{postgres_cfg.host}:{postgres_cfg.port}/{postgres_cfg.db}"
            )

        # Configurable pool settings
        pool_size = int(os.getenv("ASYNC_POOL_SIZE", "20"))
        max_overflow = int(os.getenv("ASYNC_MAX_OVERFLOW", "10"))

        _async_engine = create_async_engine(
            db_url,
            pool_size=pool_size,
            max_overflow=max_overflow,
            pool_pre_ping=True,  # Validates connections before use
            echo=postgres_cfg.echo,
        )

    return _async_engine


def get_async_session_local() -> async_sessionmaker[AsyncSession]:
    """Get or create the async session factory."""
    global _AsyncSessionLocal

    if _AsyncSessionLocal is None:
        _AsyncSessionLocal = async_sessionmaker(
            bind=get_async_engine(),
            class_=AsyncSession,
            expire_on_commit=False,
        )

    return _AsyncSessionLocal


async def get_async_db() -> AsyncGenerator[AsyncSession, None]:
    """Get an async database session."""
    SessionLocal = get_async_session_local()
    async with SessionLocal() as session:
        yield session
