# Testing Guide

**Last Updated**: 2025-12-04  
**Status**: Testing Infrastructure

## Overview

This guide covers testing strategies, test structure, and best practices for the multi-omics RAG platform.

## Test Structure

```
tests/
├── unit/
│   ├── ingestion/
│   │   ├── test_metabolomics_ingestion.py
│   │   ├── test_proteomics_ingestion.py
│   │   ├── test_transcriptomics_ingestion.py
│   │   └── test_lipidomics_ingestion.py
│   ├── parsing/
│   ├── scoring/
│   └── embedding/
├── integration/
│   ├── test_postgres_integration.py
│   ├── test_notion_sync.py
│   └── test_end_to_end.py
└── fixtures/
    ├── sample_data/
    └── mock_responses/
```

## Running Tests

### All Tests
```bash
pytest
```

### Specific Test Suite
```bash
# Unit tests only
pytest tests/unit/

# Integration tests only
pytest tests/integration/

# Specific test file
pytest tests/unit/ingestion/test_metabolomics_ingestion.py
```

### With Coverage
```bash
pytest --cov=amprenta_rag --cov-report=html
```

## Unit Tests

### Ingestion Tests

Test individual ingestion functions in isolation:

```python
# tests/unit/ingestion/test_metabolomics_ingestion.py
import pytest
from pathlib import Path
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file

def test_ingest_metabolomics_file_success(tmp_path):
    """Test successful metabolomics ingestion."""
    # Create sample file
    sample_file = tmp_path / "sample.csv"
    sample_file.write_text("Metabolite,Value\nGlucose,1.5\nLactate,2.3")
    
    # Mock Postgres and Notion
    # ... test ingestion ...
    
def test_ingest_metabolomics_file_missing_file():
    """Test error handling for missing file."""
    with pytest.raises(FileNotFoundError):
        ingest_metabolomics_file("nonexistent.csv")
```

### Parsing Tests

Test file parsing functions:

```python
# tests/unit/parsing/test_metabolite_parsing.py
from amprenta_rag.ingestion.metabolomics.file_parsing import (
    extract_metabolite_set_from_file,
    normalize_metabolite_name,
)

def test_normalize_metabolite_name():
    """Test metabolite name normalization."""
    assert normalize_metabolite_name("Glucose") == "glucose"
    assert normalize_metabolite_name("L-lactate") == "lactate"
    assert normalize_metabolite_name("alpha-Ketoglutarate") == "ketoglutarate"
```

### Scoring Tests

Test signature scoring:

```python
# tests/unit/scoring/test_signature_scoring.py
from amprenta_rag.ingestion.multi_omics_scoring import score_signature_against_dataset

def test_score_signature_perfect_match():
    """Test signature scoring with perfect match."""
    signature_features = {"feature1", "feature2", "feature3"}
    dataset_features = {"feature1", "feature2", "feature3"}
    
    score = score_signature_against_dataset(
        signature_features=signature_features,
        dataset_features=dataset_features,
    )
    
    assert score["overlap_fraction"] == 1.0
    assert score["matched_count"] == 3
```

## Integration Tests

### Postgres Integration

Test Postgres dataset creation and retrieval:

```python
# tests/integration/test_postgres_integration.py
import pytest
from amprenta_rag.database.base import get_db
from amprenta_rag.models.domain import OmicsType
from amprenta_rag.ingestion.postgres_integration import (
    create_or_update_dataset_in_postgres,
)

@pytest.fixture
def db_session():
    """Create test database session."""
    db = next(get_db())
    yield db
    db.rollback()
    db.close()

def test_create_dataset_in_postgres(db_session):
    """Test creating dataset in Postgres."""
    dataset = create_or_update_dataset_in_postgres(
        name="Test Dataset",
        omics_type=OmicsType.METABOLOMICS,
        file_paths=["test.csv"],
        description="Test description",
    )
    
    assert dataset.id is not None
    assert dataset.name == "Test Dataset"
    assert dataset.omics_type == OmicsType.METABOLOMICS
```

### End-to-End Tests

Test complete ingestion workflow:

```python
# tests/integration/test_end_to_end.py
def test_full_metabolomics_ingestion_workflow(tmp_path):
    """Test complete metabolomics ingestion."""
    # 1. Create sample file
    sample_file = tmp_path / "sample.csv"
    sample_file.write_text("Metabolite,Value\nGlucose,1.5")
    
    # 2. Ingest file
    dataset_id = ingest_metabolomics_file(
        file_path=str(sample_file),
        create_page=False,
    )
    
    # 3. Verify Postgres dataset
    db = next(get_db())
    dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    assert dataset is not None
    assert dataset.omics_type == OmicsType.METABOLOMICS
    
    # 4. Verify embedding (check Pinecone)
    # ...
```

## Test Fixtures

### Sample Data Files

Create reusable test fixtures:

```python
# tests/fixtures/sample_data.py
import pytest

@pytest.fixture
def sample_metabolomics_file(tmp_path):
    """Create sample metabolomics CSV file."""
    file_path = tmp_path / "metabolomics_sample.csv"
    file_path.write_text("""
Metabolite,Value,StdDev
Glucose,1.5,0.1
Lactate,2.3,0.2
Pyruvate,0.8,0.05
    """.strip())
    return str(file_path)

@pytest.fixture
def sample_proteomics_file(tmp_path):
    """Create sample proteomics CSV file."""
    file_path = tmp_path / "proteomics_sample.csv"
    file_path.write_text("""
Protein,Log2FC,PValue
P12345,1.2,0.001
Q67890,-0.8,0.05
    """.strip())
    return str(file_path)
```

### Mock Services

Mock external services:

```python
# tests/fixtures/mocks.py
import pytest
from unittest.mock import Mock, patch

@pytest.fixture
def mock_notion_client():
    """Mock Notion API client."""
    with patch("amprenta_rag.clients.notion_client") as mock:
        yield mock

@pytest.fixture
def mock_pinecone():
    """Mock Pinecone client."""
    with patch("amprenta_rag.clients.pinecone_client") as mock:
        yield mock
```

## Test Configuration

### pytest.ini

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = 
    -v
    --tb=short
    --strict-markers
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow tests
    requires_postgres: Tests requiring Postgres
    requires_notion: Tests requiring Notion
```

### Test Environment

Create `.env.test`:

```bash
# Test database (separate from production)
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=amprenta_rag_test
POSTGRES_USER=test_user
POSTGRES_PASSWORD=test_password

# Use Postgres in tests
USE_POSTGRES_AS_SOT=true
ENABLE_NOTION_SYNC=false

# Test API keys (safe to commit)
NOTION_API_KEY=test_key
PINECONE_API_KEY=test_key
```

## Best Practices

### 1. Isolation
- Each test should be independent
- Use fixtures for setup/teardown
- Clean up after tests

### 2. Speed
- Use mocks for external services
- Separate fast unit tests from slow integration tests
- Use markers to skip slow tests: `@pytest.mark.slow`

### 3. Coverage
- Aim for >80% code coverage
- Focus on critical paths
- Test error cases

### 4. Documentation
- Use descriptive test names
- Add docstrings to test functions
- Document test fixtures

## Continuous Integration

### GitHub Actions Example

```yaml
# .github/workflows/test.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    services:
      postgres:
        image: postgres:14
        env:
          POSTGRES_DB: amprenta_rag_test
          POSTGRES_USER: test_user
          POSTGRES_PASSWORD: test_password
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5
    
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install pytest pytest-cov
      
      - name: Run tests
        run: pytest --cov=amprenta_rag --cov-report=xml
      
      - name: Upload coverage
        uses: codecov/codecov-action@v2
```

## Test Checklist

### For Each New Feature

- [ ] Unit tests for core functions
- [ ] Integration tests for workflows
- [ ] Error handling tests
- [ ] Edge case tests
- [ ] Documentation updated

### Before Release

- [ ] All tests passing
- [ ] Coverage >80%
- [ ] No slow tests in default suite
- [ ] CI/CD pipeline passing

## Resources

- pytest documentation: https://docs.pytest.org/
- Coverage.py: https://coverage.readthedocs.io/
- Testing best practices: See `docs/DEVELOPMENT.md`

