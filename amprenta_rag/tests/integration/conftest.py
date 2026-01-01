"""Integration test fixtures with real database."""

import pytest

# P1 FIX: Auto-apply integration marker to all tests in this directory
pytestmark = pytest.mark.integration

from uuid import uuid4
from sqlalchemy.orm import Session
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.base import get_session_local, get_engine
from amprenta_rag.database.models import User, Program


@pytest.fixture(scope="module")
def integration_db():
    """Get real database session for integration tests."""
    SessionLocal = get_session_local()
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@pytest.fixture
def db_session(integration_db):
    """Transaction-scoped session with rollback for test isolation."""
    connection = get_engine().connect()
    transaction = connection.begin()
    session = Session(bind=connection)
    
    yield session
    
    session.close()
    transaction.rollback()
    connection.close()


@pytest.fixture
def test_user(db_session):
    """Create a real test user in database."""
    user = User(
        id=uuid4(),
        username=f"testuser_{uuid4().hex[:8]}",
        email=f"test_{uuid4().hex[:8]}@test.com",
        hashed_password="fakehash",
        is_active=True,
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def test_program(db_session):
    """Create a real test program in database."""
    program = Program(
        id=uuid4(),
        name=f"Test Program {uuid4().hex[:8]}",
        description="Integration test program",
    )
    db_session.add(program)
    db_session.commit()
    db_session.refresh(program)
    return program


@pytest.fixture
def integration_client(db_session, test_user):
    """TestClient with real DB session and authenticated user."""
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    
    def override_get_db():
        yield db_session
    
    def override_get_user():
        return test_user
    
    app.dependency_overrides[get_db] = override_get_db
    app.dependency_overrides[get_current_user] = override_get_user
    
    client = TestClient(app)
    yield client
    
    app.dependency_overrides.clear()
