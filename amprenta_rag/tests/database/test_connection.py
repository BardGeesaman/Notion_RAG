"""
Tests for Postgres database connection.

These tests verify:
- Database connection works
- Configuration is correct
- Connection pooling functions
"""

import pytest
from sqlalchemy import text

from amprenta_rag.config import get_config
from amprenta_rag.database.base import get_engine, get_session_local


class TestDatabaseConnection:
    """Test database connection functionality."""
    
    def test_config_loaded(self):
        """Test that Postgres configuration can be loaded."""
        cfg = get_config()
        postgres_cfg = cfg.postgres
        
        assert postgres_cfg is not None
        assert postgres_cfg.host is not None
        assert postgres_cfg.db is not None
        assert postgres_cfg.user is not None
    
    def test_engine_creation(self):
        """Test that database engine can be created."""
        engine = get_engine()
        assert engine is not None
        assert str(engine.url).startswith("postgresql")
    
    def test_connection(self):
        """Test that we can connect to the database."""
        engine = get_engine()
        
        with engine.connect() as conn:
            result = conn.execute(text("SELECT 1"))
            assert result.scalar() == 1
    
    def test_version_query(self):
        """Test that we can query PostgreSQL version."""
        engine = get_engine()
        
        with engine.connect() as conn:
            result = conn.execute(text("SELECT version()"))
            version = result.scalar()
            assert version is not None
            assert "PostgreSQL" in version
    
    def test_session_creation(self):
        """Test that database sessions can be created."""
        SessionLocal = get_session_local()
        assert SessionLocal is not None
        
        db = SessionLocal()
        try:
            # Test a simple query
            result = db.execute(text("SELECT 1"))
            assert result.scalar() == 1
        finally:
            db.close()


class TestDatabaseSchema:
    """Test database schema exists."""
    
    def test_tables_exist(self):
        """Test that expected tables exist in the database."""
        engine = get_engine()
        
        expected_tables = [
            "programs",
            "experiments",
            "datasets",
            "features",
            "signatures",
            "signature_components",
        ]
        
        with engine.connect() as conn:
            result = conn.execute(
                text("""
                    SELECT table_name 
                    FROM information_schema.tables 
                    WHERE table_schema = 'public'
                    ORDER BY table_name
                """)
            )
            existing_tables = {row[0] for row in result}
        
        for table in expected_tables:
            assert table in existing_tables, f"Table '{table}' not found in database"
    
    def test_alembic_version_table(self):
        """Test that Alembic version table exists (migrations applied)."""
        engine = get_engine()
        
        with engine.connect() as conn:
            result = conn.execute(
                text("""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables 
                        WHERE table_schema = 'public' 
                        AND table_name = 'alembic_version'
                    )
                """)
            )
            exists = result.scalar()
        
        assert exists, "Alembic version table not found - migrations may not be applied"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

