"""Pytest configuration and fixtures for all tests.

Provides common fixtures and setup for database and API tests.
"""

# MUST be first - configure test environment before ANY app imports
from amprenta_rag.test_config import configure_test_environment
configure_test_environment()

import os
import pytest
from pathlib import Path


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "postgres: marks tests that require Postgres database"
    )
    config.addinivalue_line(
        "markers", "api: marks tests that require FastAPI server"
    )
    config.addinivalue_line(
        "markers", "slow: marks tests as slow running"
    )


@pytest.fixture(scope="session")
def requires_postgres():
    """Check if Postgres is available, skip test if not."""
    try:
        from amprenta_rag.config import get_config
        cfg = get_config()
        postgres_cfg = cfg.postgres

        # Check if Postgres config is present
        if not postgres_cfg.user or not postgres_cfg.db:
            pytest.skip("Postgres not configured (set POSTGRES_USER, POSTGRES_DB)")

        # Try to connect
        from amprenta_rag.database.base import get_engine
        engine = get_engine()
        with engine.connect() as conn:
            conn.execute("SELECT 1")

        return True
    except Exception as e:
        pytest.skip(f"Postgres not available: {e}")


@pytest.fixture(scope="session")
def requires_fastapi():
    """Check if FastAPI is available, skip test if not."""
    try:
        return True
    except ImportError:
        pytest.skip("FastAPI not installed")


@pytest.fixture(scope="session")
def requires_psycopg2():
    """Check if psycopg2 is available, skip test if not."""
    try:
        return True
    except ImportError:
        pytest.skip("psycopg2 not installed (required for Postgres)")


@pytest.fixture(scope="session", autouse=True)
def ensure_migrations_applied():
    """
    Ensure Alembic migrations are applied before running tests.
    
    This fixture runs automatically for all tests (autouse=True) at session scope.
    It automatically applies any pending migrations to ensure the test database
    schema is up to date.
    """
    print("\n" + "="*70)
    print("RUNNING MIGRATION FIXTURE")
    print("="*70)
    
    try:
        from alembic import command
        from alembic.config import Config
        from amprenta_rag.database.base import get_engine
        from sqlalchemy import text
        import os
        import sys
        
        # Check if database is available
        print("Step 1: Checking database availability...")
        try:
            engine = get_engine()
            with engine.connect() as conn:
                result = conn.execute(text("SELECT 1"))
                print(f"  ✓ Database connected: {engine.url}")
        except Exception as e:
            print(f"  ✗ Database not available: {e}")
            print("  → Tests that need DB will fail naturally")
            return
        
        # Set up Alembic configuration
        print("\nStep 2: Locating alembic.ini...")
        alembic_ini_path = Path(__file__).parent.parent.parent / "alembic.ini"
        print(f"  Trying: {alembic_ini_path}")
        
        if not alembic_ini_path.exists():
            # Try from current working directory
            alembic_ini_path = Path.cwd() / "alembic.ini"
            print(f"  Trying: {alembic_ini_path}")
        
        if not alembic_ini_path.exists():
            raise FileNotFoundError(f"alembic.ini not found. Searched: {alembic_ini_path}")
        
        print(f"  ✓ Found: {alembic_ini_path}")
        
        # Set up Alembic configuration
        print("\nStep 3: Configuring Alembic...")
        alembic_cfg = Config(str(alembic_ini_path))
        
        # Get script location from ini file or default to "alembic"
        script_location = alembic_cfg.get_main_option("script_location", "alembic")
        print(f"  Script location: {script_location}")
        
        # Override sqlalchemy.url to use the same engine as tests
        alembic_cfg.set_main_option("sqlalchemy.url", str(engine.url))
        print(f"  Database URL: {engine.url}")
        
        # Run migrations to head
        print("\nStep 4: Running migrations (alembic upgrade head)...")
        
        # Capture alembic output
        from io import StringIO
        import logging
        
        # Configure alembic logging to see what's happening
        alembic_logger = logging.getLogger('alembic')
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        alembic_logger.addHandler(handler)
        alembic_logger.setLevel(logging.INFO)
        
        command.upgrade(alembic_cfg, "head")
        
        print("  ✓ Migrations completed successfully")
        print("="*70 + "\n")
        
    except Exception as e:
        print(f"\n  ✗ MIGRATION FAILED: {e}")
        print(f"  Type: {type(e).__name__}")
        import traceback
        print("  Traceback:")
        traceback.print_exc()
        print("="*70 + "\n")
        
        # RAISE the error instead of warning
        # This will cause tests to fail immediately with clear error
        raise RuntimeError(f"Failed to apply database migrations: {e}") from e

