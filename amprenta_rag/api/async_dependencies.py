"""FastAPI async dependencies for database sessions."""

from __future__ import annotations

from typing import AsyncGenerator

from sqlalchemy.ext.asyncio import AsyncSession

from amprenta_rag.database.async_base import get_async_db


async def get_async_database_session() -> AsyncGenerator[AsyncSession, None]:
    """Dependency to get an async database session."""
    async for session in get_async_db():
        yield session
