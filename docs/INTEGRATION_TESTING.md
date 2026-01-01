# Integration Testing Guide

Last Updated: 2025-01-01

## Overview

This guide documents the integration test infrastructure and methodology for converting mock-heavy API tests to real database tests. The goal is to improve test reliability and catch integration issues while maintaining fast execution times through strategic mocking.

## Infrastructure

### Fixtures (`amprenta_rag/tests/integration/conftest.py`)

The integration test infrastructure provides comprehensive fixtures for real database testing:

#### Core Database Fixtures
- `integration_db` - Module-scoped real database session
- `db_session` - Transaction-scoped session with automatic rollback for test isolation
- `integration_client` - FastAPI TestClient with dependency overrides

#### Entity Fixtures
- `test_user` - Real user entity with unique credentials
- `test_program` - Real program entity linked to test user
- `test_compound` - Real compound entity for chemistry tests
- `test_dataset` - Real dataset entity for omics tests
- `test_experiment` - Real experiment entity with program linkage
- `test_signature` - Real signature entity for analysis tests

#### Authentication Fixtures
- `admin_user` - Admin user with elevated privileges
- `admin_client` - TestClient with admin authentication pre-configured

#### Performance Fixtures
- `timed_request` - Performance-tracked API calls with automatic benchmark recording
- `benchmark_tracker` - Performance metric collection and threshold validation

### Benchmarks (`amprenta_rag/tests/integration/benchmarks.py`)

Performance tracking infrastructure with configurable thresholds:

#### BenchmarkTracker Class
- Tracks operation timing with method-specific thresholds
- Generates JSON reports for CI integration
- Provides pass/fail status based on performance gates

#### Default Thresholds
- `GET`: 150ms - Simple read operations
- `POST`: 200ms - Create operations
- `PUT/PATCH`: 200ms - Update operations  
- `DELETE`: 100ms - Delete operations
- `batch`: 500ms - Batch operations
- `compute`: 3000ms - ML/compute operations (ADMET, QSAR)

## Conversion Methodology

### Step 1: Create Real Entities

Replace mocked database objects with real entities:

```python
# Before (mock-heavy)
@patch("amprenta_rag.database.session.db_session")
def test_create_job(mock_db):
    mock_job = MagicMock()
    mock_db.return_value.__enter__.return_value.query.return_value.first.return_value = mock_job
    # ... test logic

# After (integration)
@pytest.mark.integration
def test_create_job_integration(integration_client, db_session, timed_request):
    job = Job(id=uuid4(), type="export", status="pending", created_at=datetime.utcnow())
    db_session.add(job)
    db_session.commit()
    
    response, benchmark = timed_request("POST", f"/api/v1/jobs/{job.id}/trigger", "test_create_job")
    assert response.status_code == 200
    assert benchmark.passed
```

### Step 2: Mock Only External Services

Keep real database interactions while mocking external dependencies:

```python
@pytest.mark.integration
@patch("amprenta_rag.api.routers.jobs.celery_app")  # Mock Celery
@patch("amprenta_rag.services.backup.s3_client")    # Mock S3
def test_backup_operation(mock_s3, mock_celery, integration_client, db_session):
    # Real database entities
    backup = BackupRecord(id=uuid4(), status="pending")
    db_session.add(backup)
    db_session.commit()
    
    # Mock external services
    mock_celery.send_task.return_value.id = "task-123"
    mock_s3.upload_file.return_value = True
    
    # Test with real DB + mocked externals
    response = integration_client.post(f"/api/v1/backup/{backup.id}/execute")
    assert response.status_code == 202
```

### Step 3: Use Performance Assertions

Integrate performance tracking into all tests:

```python
@pytest.mark.integration
def test_list_compounds(integration_client, timed_request):
    response, benchmark = timed_request(
        "GET", 
        "/api/v1/compounds", 
        "test_list_compounds",
        params={"limit": 50}
    )
    
    assert response.status_code == 200
    assert benchmark.passed, f"Took {benchmark.elapsed_ms}ms (threshold: {benchmark.threshold_ms}ms)"
    assert len(response.json()) <= 50
```

### Step 4: Archive Original

Move original mock-heavy tests to preserve history:

```bash
# Archive pattern
mv amprenta_rag/tests/api/test_jobs_api.py amprenta_rag/tests/api/_archived/
```

## What to Mock vs Real

| Category | Approach | Rationale |
|----------|----------|-----------|
| **Database CRUD** | Real | Core business logic, catch FK violations |
| **FK relationships** | Real | Verify referential integrity |
| **Status transitions** | Real | Critical state machine logic |
| **User authentication** | Real | Security-critical flows |
| **External APIs** | Mock | Unreliable, slow, rate-limited |
| **File storage (S3)** | Mock | Expensive, slow, requires credentials |
| **Celery execution** | Mock | Background tasks, timing issues |
| **Email sending** | Mock | External service, no test value |
| **ML model inference** | Mock | Slow, GPU-dependent |
| **File system ops** | `tempfile` | Real but isolated |

## Running Integration Tests

### Basic Execution

```bash
# Run all integration tests
pytest -m integration -v

# Run specific domain
pytest amprenta_rag/tests/integration/api/test_jobs_integration.py -v

# Run with performance tracking
pytest -m integration --tb=short
```

### Performance Monitoring

```bash
# Generate benchmark report
pytest -m integration --json-report --json-report-file=benchmark_results.json

# Check performance gates (CI)
python scripts/check_benchmarks.py benchmark_results.json

# Allow higher failure rate for development
python scripts/check_benchmarks.py benchmark_results.json --max-fail 0.1
```

### CI Integration

Add to GitHub Actions workflow:

```yaml
- name: Run integration tests with benchmarks
  run: |
    pytest -m integration --json-report --json-report-file=benchmark_results.json
    python scripts/check_benchmarks.py benchmark_results.json
```

## Troubleshooting

### Common Issues

#### Foreign Key Violations
```python
# Problem: Random UUID for FK
log_activity(actor_id=uuid4())  # FK violation!

# Solution: Use None or real entity
log_activity(actor_id=None)  # System action
# OR
user = User(id=uuid4(), username=f"test_{uuid4().hex[:8]}")
db_session.add(user)
db_session.commit()
log_activity(actor_id=user.id)  # Valid FK
```

#### Unique Constraint Violations
```python
# Problem: Hardcoded values
user = User(username="testuser", email="test@test.com")  # Collision!

# Solution: UUID-based unique values
user = User(
    username=f"testuser_{uuid4().hex[:8]}", 
    email=f"test_{uuid4().hex[:8]}@test.com"
)
```

#### DetachedInstanceError
```python
# Problem: Accessing object outside session
def get_user():
    with db_session() as db:
        user = db.query(User).first()
        return user  # Detached!

# Solution: Expunge before returning
def get_user():
    with db_session() as db:
        user = db.query(User).first()
        if user:
            db.expunge(user)  # Critical!
        return user
```

### Performance Issues

#### Slow Test Identification
```bash
# Find slowest tests
pytest -m integration --durations=10

# Check benchmark report
python scripts/check_benchmarks.py --max-fail 0.0  # Show all failures
```

#### Threshold Tuning
Adjust thresholds in `BenchmarkTracker.THRESHOLDS` based on actual performance:

```python
THRESHOLDS = {
    "GET": 200,      # Increased for complex queries
    "POST": 300,     # Increased for validation-heavy creates
    "compute": 5000, # Increased for ML operations
}
```

## Metrics & Success Criteria

### Conversion Metrics
- **Files Converted**: 6 (Phase 1-2)
- **Integration Tests**: 81 total
- **Mock Reduction**: 87% average
- **Domains Covered**: Jobs, Imaging, Backup, Sync, Papers, Collaboration

### Quality Gates
- **Regression Tests**: 0 (all tests pass)
- **Performance Failures**: <5% of tests exceed thresholds
- **Linting Errors**: 0 
- **Coverage Maintained**: ✅

### Infrastructure Readiness
- **Fixtures Available**: 12 comprehensive fixtures
- **Benchmark Framework**: Operational
- **CI Integration**: GitHub Actions ready
- **Team Enablement**: Documentation complete

## Future Expansion

### Phase 2 Domains (Ready for Conversion)
- **Medium-Risk**: Papers citations, Teams, IP, Protocols, Extraction, Screening, Annotations (8 files)
- **Low-Risk**: Comments, Audit, Portfolio, Export, Activity, Reviews, Sharing (8 files)

### Phase 3 Domains (Planned)
- **ML/Analysis**: ADMET, QSAR, Signatures, Features (15 files)
- **Chemistry/HTS**: Compounds, Campaigns, Plates, Wells (10 files)  
- **Core CRUD**: Programs, Datasets, Experiments, Users (15 files)

### Methodology Scaling
The proven conversion methodology can be applied to any API test file:

1. **Analyze**: Count mocks, identify external dependencies
2. **Convert**: Real DB entities, mock externals, add performance tracking
3. **Verify**: All tests pass, performance within thresholds
4. **Archive**: Move original to `_archived/`

Each conversion follows the same pattern, making large-scale adoption straightforward.

## References

- **Baseline Report**: `reports/integration_test_baseline.json`
- **CI Script**: `scripts/check_benchmarks.py`
- **Fixtures**: `amprenta_rag/tests/integration/conftest.py`
- **Benchmarks**: `amprenta_rag/tests/integration/benchmarks.py`
- **Examples**: `amprenta_rag/tests/integration/api/test_*_integration.py`
