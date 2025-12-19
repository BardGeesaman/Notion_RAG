"""
Pytest configuration and fixtures for all tests.

Provides common fixtures and setup for database and API tests.
"""

import pytest


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

