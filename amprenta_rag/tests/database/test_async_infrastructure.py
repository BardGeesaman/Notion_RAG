"""Tests for async database infrastructure."""

from __future__ import annotations

from unittest.mock import MagicMock, patch
from uuid import uuid4

import pytest
from sqlalchemy.ext.asyncio import AsyncEngine, AsyncSession, async_sessionmaker


class TestAsyncEngine:
    """Test async engine creation and configuration."""

    @pytest.mark.asyncio
    async def test_async_engine_creation(self):
        """Verify async engine initializes with correct settings."""
        with patch('amprenta_rag.database.async_base.get_config') as mock_config:
            with patch('amprenta_rag.database.async_base.create_async_engine') as mock_create_engine:
                # Mock config
                mock_postgres = MagicMock()
                mock_postgres.url = "postgresql://user:pass@localhost:5432/testdb"
                mock_postgres.echo = False
                mock_config.return_value.postgres = mock_postgres
                
                # Mock engine
                mock_engine = MagicMock(spec=AsyncEngine)
                mock_create_engine.return_value = mock_engine
                
                from amprenta_rag.database.async_base import get_async_engine
                
                # Clear any existing engine
                import amprenta_rag.database.async_base as async_base_module
                async_base_module._async_engine = None
                
                engine = get_async_engine()
                
                # Verify engine creation was called with correct parameters
                mock_create_engine.assert_called_once()
                call_args, call_kwargs = mock_create_engine.call_args
                
                assert call_args[0] == "postgresql+asyncpg://user:pass@localhost:5432/testdb"
                assert call_kwargs["pool_size"] == 20  # Default value
                assert call_kwargs["max_overflow"] == 10  # Default value
                assert call_kwargs["pool_pre_ping"] == True
                assert call_kwargs["echo"] == False
                
                assert engine == mock_engine

    @pytest.mark.asyncio
    async def test_async_engine_with_env_vars(self):
        """Verify async engine respects environment variables."""
        with patch('amprenta_rag.database.async_base.get_config') as mock_config:
            with patch('amprenta_rag.database.async_base.create_async_engine') as mock_create_engine:
                with patch('amprenta_rag.database.async_base.os.getenv') as mock_getenv:
                    # Mock config
                    mock_postgres = MagicMock()
                    mock_postgres.url = None
                    mock_postgres.user = "testuser"
                    mock_postgres.password = "testpass"
                    mock_postgres.host = "testhost"
                    mock_postgres.port = 5432
                    mock_postgres.db = "testdb"
                    mock_postgres.echo = True
                    mock_config.return_value.postgres = mock_postgres
                    
                    # Mock environment variables
                    def getenv_side_effect(key, default):
                        env_vars = {
                            "ASYNC_POOL_SIZE": "30",
                            "ASYNC_MAX_OVERFLOW": "15"
                        }
                        return env_vars.get(key, default)
                    
                    mock_getenv.side_effect = getenv_side_effect
                    
                    # Mock engine
                    mock_engine = MagicMock(spec=AsyncEngine)
                    mock_create_engine.return_value = mock_engine
                    
                    from amprenta_rag.database.async_base import get_async_engine
                    
                    # Clear any existing engine
                    import amprenta_rag.database.async_base as async_base_module
                    async_base_module._async_engine = None
                    
                    engine = get_async_engine()
                    
                    # Verify engine creation was called with environment variables
                    mock_create_engine.assert_called_once()
                    call_args, call_kwargs = mock_create_engine.call_args
                    
                    assert call_args[0] == "postgresql+asyncpg://testuser:testpass@testhost:5432/testdb"
                    assert call_kwargs["pool_size"] == 30  # From env var
                    assert call_kwargs["max_overflow"] == 15  # From env var
                    assert call_kwargs["pool_pre_ping"] == True
                    assert call_kwargs["echo"] == True


class TestAsyncSessionFactory:
    """Test async session factory creation."""

    @pytest.mark.asyncio
    async def test_async_session_factory(self):
        """Verify async session factory creates valid sessions."""
        with patch('amprenta_rag.database.async_base.get_async_engine') as mock_get_engine:
            with patch('amprenta_rag.database.async_base.async_sessionmaker') as mock_sessionmaker:
                # Mock engine
                mock_engine = MagicMock(spec=AsyncEngine)
                mock_get_engine.return_value = mock_engine
                
                # Mock session factory
                mock_factory = MagicMock(spec=async_sessionmaker)
                mock_sessionmaker.return_value = mock_factory
                
                from amprenta_rag.database.async_base import get_async_session_local
                
                # Clear any existing session factory
                import amprenta_rag.database.async_base as async_base_module
                async_base_module._AsyncSessionLocal = None
                
                factory = get_async_session_local()
                
                # Verify session factory creation
                mock_sessionmaker.assert_called_once_with(
                    bind=mock_engine,
                    class_=AsyncSession,
                    expire_on_commit=False,
                )
                
                assert factory == mock_factory


class TestAsyncContextManager:
    """Test async database session context manager."""

    @pytest.mark.asyncio
    async def test_async_context_manager_commit(self):
        """Verify context manager commits on success."""
        with patch('amprenta_rag.database.async_session.get_async_session_local') as mock_get_factory:
            # Mock session factory and session with async methods
            mock_session = MagicMock(spec=AsyncSession)
            
            # Make commit and rollback async
            async def mock_commit():
                pass
            async def mock_rollback():
                pass
            
            mock_session.commit = MagicMock(side_effect=mock_commit)
            mock_session.rollback = MagicMock(side_effect=mock_rollback)
            
            mock_factory = MagicMock()
            mock_factory.return_value.__aenter__.return_value = mock_session
            mock_factory.return_value.__aexit__.return_value = None
            mock_get_factory.return_value = mock_factory
            
            from amprenta_rag.database.async_session import async_db_session
            
            # Use the context manager successfully
            async with async_db_session() as session:
                assert session == mock_session
                # Simulate some database work
                pass
            
            # Verify commit was called, rollback was not
            mock_session.commit.assert_called_once()
            mock_session.rollback.assert_not_called()

    @pytest.mark.asyncio
    async def test_async_context_manager_rollback(self):
        """Verify context manager rolls back on exception."""
        with patch('amprenta_rag.database.async_session.get_async_session_local') as mock_get_factory:
            # Mock session factory and session with async methods
            mock_session = MagicMock(spec=AsyncSession)
            
            # Make commit and rollback async
            async def mock_commit():
                pass
            async def mock_rollback():
                pass
            
            mock_session.commit = MagicMock(side_effect=mock_commit)
            mock_session.rollback = MagicMock(side_effect=mock_rollback)
            
            mock_factory = MagicMock()
            mock_factory.return_value.__aenter__.return_value = mock_session
            mock_factory.return_value.__aexit__.return_value = None
            mock_get_factory.return_value = mock_factory
            
            from amprenta_rag.database.async_session import async_db_session
            
            # Use the context manager with an exception
            test_exception = ValueError("Test error")
            with pytest.raises(ValueError, match="Test error"):
                async with async_db_session() as session:
                    assert session == mock_session
                    raise test_exception
            
            # Verify rollback was called, commit was not
            mock_session.rollback.assert_called_once()
            mock_session.commit.assert_not_called()


class TestAsyncDependencyInjection:
    """Test FastAPI async dependency injection."""

    @pytest.mark.asyncio
    async def test_async_dependency_injection(self):
        """Verify FastAPI dependency injection works."""
        with patch('amprenta_rag.api.async_dependencies.get_async_db') as mock_get_db:
            # Mock session
            mock_session = MagicMock(spec=AsyncSession)
            
            # Mock async generator
            async def mock_async_generator():
                yield mock_session
            
            mock_get_db.return_value = mock_async_generator()
            
            from amprenta_rag.api.async_dependencies import get_async_database_session
            
            # Test the dependency
            session_generator = get_async_database_session()
            session = await session_generator.__anext__()
            
            assert session == mock_session
            
            # Verify the generator is exhausted
            with pytest.raises(StopAsyncIteration):
                await session_generator.__anext__()
