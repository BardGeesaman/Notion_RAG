"""Integration test fixtures with real database."""

import pytest
import time
import logging
from contextlib import contextmanager

# P1 FIX: Auto-apply integration marker to all tests in this directory
pytestmark = pytest.mark.integration

from uuid import uuid4
from sqlalchemy.orm import Session
from fastapi.testclient import TestClient

from amprenta_rag.api.main import app
from amprenta_rag.database.base import get_session_local, get_engine
from amprenta_rag.database.models import (
    Program, Compound, Dataset, Experiment, Signature, Feature
)
from amprenta_rag.models.auth import User

logger = logging.getLogger(__name__)


# Performance tracking
@contextmanager
def track_performance(operation: str, threshold_ms: float = 200):
    """Track and warn on slow operations."""
    start = time.perf_counter()
    yield
    elapsed_ms = (time.perf_counter() - start) * 1000
    if elapsed_ms > threshold_ms:
        logger.warning(f"SLOW: {operation} took {elapsed_ms:.1f}ms (threshold: {threshold_ms}ms)")


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
        password_hash="fakehash",
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


# Domain fixtures
@pytest.fixture
def test_compound(db_session, test_program):
    """Create real compound with program link."""
    compound = Compound(
        id=uuid4(),
        smiles="CCO",  # ethanol
        name=f"Test Compound {uuid4().hex[:8]}",
        program_id=test_program.id,
    )
    db_session.add(compound)
    db_session.commit()
    db_session.refresh(compound)
    return compound


@pytest.fixture
def test_dataset(db_session, test_program):
    """Create real dataset with program link."""
    dataset = Dataset(
        id=uuid4(),
        name=f"Test Dataset {uuid4().hex[:8]}",
        program_id=test_program.id,
        omics_type="transcriptomics",
    )
    db_session.add(dataset)
    db_session.commit()
    db_session.refresh(dataset)
    return dataset


@pytest.fixture
def test_experiment(db_session, test_program, test_user):
    """Create real experiment."""
    experiment = Experiment(
        id=uuid4(),
        name=f"Test Experiment {uuid4().hex[:8]}",
        program_id=test_program.id,
        created_by_id=test_user.id,
    )
    db_session.add(experiment)
    db_session.commit()
    db_session.refresh(experiment)
    return experiment


@pytest.fixture
def test_signature(db_session, test_dataset):
    """Create real signature."""
    signature = Signature(
        id=uuid4(),
        name=f"Test Signature {uuid4().hex[:8]}",
        dataset_id=test_dataset.id,
    )
    db_session.add(signature)
    db_session.commit()
    db_session.refresh(signature)
    return signature


@pytest.fixture
def admin_user(db_session):
    """Create admin user for admin endpoint tests."""
    user = User(
        id=uuid4(),
        username=f"admin_{uuid4().hex[:8]}",
        email=f"admin_{uuid4().hex[:8]}@test.com",
        password_hash="fakehash",
        is_active=True,
        role="admin",
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture
def admin_client(db_session, admin_user):
    """TestClient with admin authentication."""
    from amprenta_rag.database.base import get_db
    from amprenta_rag.api.dependencies import get_current_user
    from amprenta_rag.api.main import app
    from fastapi.testclient import TestClient
    
    def override_get_db():
        yield db_session
    
    def override_get_user():
        return admin_user
    
    app.dependency_overrides[get_db] = override_get_db
    app.dependency_overrides[get_current_user] = override_get_user
    
    client = TestClient(app)
    yield client
    
    app.dependency_overrides.clear()


# Benchmark fixtures
from amprenta_rag.tests.integration.benchmarks import get_tracker, BenchmarkTracker


@pytest.fixture
def benchmark_tracker() -> BenchmarkTracker:
    """Get benchmark tracker for performance assertions."""
    return get_tracker()


@pytest.fixture
def timed_request(integration_client, benchmark_tracker):
    """Helper for making timed API requests."""
    def _make_request(method: str, endpoint: str, test_name: str, **kwargs):
        start = time.perf_counter()
        
        if method.upper() == "GET":
            response = integration_client.get(endpoint, **kwargs)
        elif method.upper() == "POST":
            response = integration_client.post(endpoint, **kwargs)
        elif method.upper() == "PUT":
            response = integration_client.put(endpoint, **kwargs)
        elif method.upper() == "PATCH":
            response = integration_client.patch(endpoint, **kwargs)
        elif method.upper() == "DELETE":
            response = integration_client.delete(endpoint, **kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        elapsed_ms = (time.perf_counter() - start) * 1000
        result = benchmark_tracker.record(test_name, endpoint, method, elapsed_ms)
        
        return response, result
    
    return _make_request
