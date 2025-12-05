# TIER 3 Testing & Validation Guide

**Status**: Complete  
**Last Updated**: 2025-01-XX

## Overview

This guide describes the comprehensive test suite and validation scripts created for TIER 3 infrastructure (Postgres + FastAPI). These tools verify that the database, API, and integration layers are working correctly.

## Test Suite Components

### 1. Validation Scripts

#### `scripts/validate_postgres_setup.py`
Validates Postgres database setup and connection.

**Checks**:
- Postgres configuration
- Database connectivity
- Schema tables exist
- Migrations are applied
- SQLAlchemy models are queryable

**Usage**:
```bash
python scripts/validate_postgres_setup.py
```

**Example Output**:
```
================================================================================
POSTGRES DATABASE VALIDATION
================================================================================

1. Checking Postgres Configuration...
  ✅ Configuration looks good

2. Testing Database Connection...
  ✅ Connected successfully
  PostgreSQL version: PostgreSQL 15.0...

3. Checking Schema Tables...
  ✅ All 10 expected tables exist

4. Checking Migration Status...
  ✅ Migrations applied
  Current revision: abc123def456

5. Testing SQLAlchemy Models...
  ✅ All models are queryable

✅ All checks passed! Postgres is ready to use.
```

#### `scripts/test_api_endpoints.py`
Tests FastAPI endpoints to verify API functionality.

**Checks**:
- API server is running
- All endpoints are accessible
- CRUD operations work
- API documentation is available

**Usage**:
```bash
# Make sure API server is running first
python scripts/run_api_server.py

# In another terminal
python scripts/test_api_endpoints.py
```

**Options**:
```bash
python scripts/test_api_endpoints.py --url http://localhost:8000  # Custom URL
python scripts/test_api_endpoints.py --wait 5  # Wait 5 seconds for server
```

#### `scripts/validate_tier3_complete.py`
Comprehensive validation script that runs all checks.

**Usage**:
```bash
# Run all validations
python scripts/validate_tier3_complete.py

# Skip API tests (if server not running)
python scripts/validate_tier3_complete.py --skip-api

# Skip pytest tests
python scripts/validate_tier3_complete.py --skip-tests
```

### 2. Unit Tests

#### Database Tests (`amprenta_rag/tests/database/`)

**`test_connection.py`**: Tests Postgres connection
- Configuration loading
- Engine creation
- Connection establishment
- Version queries
- Session creation

**`test_models.py`**: Tests SQLAlchemy models
- Model imports
- Model queries
- Model creation (with rollback)

**Run**:
```bash
pytest amprenta_rag/tests/database/ -v
```

#### API Tests (`amprenta_rag/tests/api/`)

**`test_endpoints.py`**: Tests FastAPI endpoints
- Root endpoints
- Programs API
- Experiments API
- Datasets API
- Features API
- Signatures API
- API documentation

**Run**:
```bash
pytest amprenta_rag/tests/api/ -v
```

### 3. Integration Tests

#### Postgres + API Integration (`amprenta_rag/tests/integration/`)

**`test_postgres_api.py`**: Tests end-to-end workflows
- Create via API → Verify in Postgres
- Read from Postgres → Verify via API
- Update via API → Verify in Postgres
- Delete via API → Verify removed from Postgres

**Run**:
```bash
pytest amprenta_rag/tests/integration/ -v
```

## Quick Start

### 1. Validate Postgres Setup

```bash
# Check Postgres connection and schema
python scripts/validate_postgres_setup.py
```

### 2. Test API Endpoints

```bash
# Start API server (in one terminal)
python scripts/run_api_server.py

# Test endpoints (in another terminal)
python scripts/test_api_endpoints.py
```

### 3. Run All Tests

```bash
# Run all pytest tests
pytest amprenta_rag/tests/ -v

# Or run comprehensive validation
python scripts/validate_tier3_complete.py
```

## Test Coverage

### ✅ What's Tested

1. **Postgres Connection**
   - Configuration validation
   - Database connectivity
   - Connection pooling

2. **Database Schema**
   - Table existence
   - Migration status
   - Model queryability

3. **SQLAlchemy Models**
   - Model imports
   - CRUD operations
   - Relationships

4. **FastAPI Endpoints**
   - All resource endpoints
   - Error handling
   - Response schemas

5. **Integration**
   - API → Postgres workflows
   - Data persistence
   - End-to-end operations

### ⚠️ What's Not Tested (Yet)

- Load/stress testing
- Concurrent access
- Migration rollbacks
- Complex relationships
- Edge cases and error scenarios

## Troubleshooting

### Postgres Connection Fails

**Error**: `Failed to connect to database`

**Solutions**:
1. Check Postgres is running: `pg_isready`
2. Verify connection settings in `.env`
3. Check database exists: `createdb amprenta_rag`
4. Verify user permissions

### API Tests Fail

**Error**: `Connection refused`

**Solutions**:
1. Start API server: `python scripts/run_api_server.py`
2. Check server is accessible: `curl http://localhost:8000/health`
3. Use `--wait` flag: `python scripts/test_api_endpoints.py --wait 5`

### Migration Tests Fail

**Error**: `Table 'X' not found`

**Solutions**:
1. Apply migrations: `python scripts/migrate_database.py`
2. Check migration status: `alembic current`
3. Verify schema: `python scripts/validate_postgres_setup.py`

## Continuous Integration

### GitHub Actions Example

```yaml
name: TIER 3 Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    services:
      postgres:
        image: postgres:15
        env:
          POSTGRES_DB: amprenta_rag
          POSTGRES_USER: postgres
          POSTGRES_PASSWORD: postgres
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5
    
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install pytest pytest-cov
      
      - name: Run database migrations
        env:
          POSTGRES_HOST: localhost
          POSTGRES_PORT: 5432
          POSTGRES_DB: amprenta_rag
          POSTGRES_USER: postgres
          POSTGRES_PASSWORD: postgres
        run: python scripts/migrate_database.py
      
      - name: Run tests
        env:
          POSTGRES_HOST: localhost
          POSTGRES_PORT: 5432
          POSTGRES_DB: amprenta_rag
          POSTGRES_USER: postgres
          POSTGRES_PASSWORD: postgres
        run: |
          pytest amprenta_rag/tests/ -v --cov=amprenta_rag
```

## Next Steps

1. **Expand Test Coverage**
   - Add more edge case tests
   - Test error scenarios
   - Add performance benchmarks

2. **Add Load Testing**
   - Concurrent request handling
   - Database connection pooling
   - API response times

3. **Automate Testing**
   - CI/CD integration
   - Pre-commit hooks
   - Nightly test runs

4. **Monitoring**
   - Test result tracking
   - Coverage reports
   - Performance metrics

## Summary

✅ **Complete Test Suite Created**:
- Validation scripts for quick checks
- Unit tests for database and API
- Integration tests for end-to-end workflows
- Comprehensive validation script

✅ **Easy to Use**:
- Simple command-line tools
- Clear error messages
- Comprehensive documentation

✅ **Production Ready**:
- Verifies infrastructure works
- Catches issues early
- Provides confidence in deployment

