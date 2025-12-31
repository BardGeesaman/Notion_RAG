"""Async database session context manager."""

from __future__ import annotations

from contextlib import asynccontextmanager
from typing import AsyncGenerator

from sqlalchemy.ext.asyncio import AsyncSession

from amprenta_rag.database.async_base import get_async_session_local


@asynccontextmanager
async def async_db_session() -> AsyncGenerator[AsyncSession, None]:
    """Async context manager for database sessions with automatic cleanup."""
    SessionLocal = get_async_session_local()
    async with SessionLocal() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
