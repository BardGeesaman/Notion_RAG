# Tests

This directory contains all tests for the multi-omics RAG platform.

## Structure

```
tests/
├── unit/              # Unit tests (fast, isolated)
│   ├── ingestion/     # Ingestion pipeline tests
│   ├── parsing/       # File parsing tests
│   ├── scoring/       # Signature scoring tests
│   └── embedding/     # RAG embedding tests
├── integration/       # Integration tests (require external services)
│   ├── postgres/      # Postgres integration tests
│   ├── notion/        # Notion sync tests
│   └── end_to_end/    # Full workflow tests
└── fixtures/          # Test data and mocks
    ├── sample_data/   # Sample CSV/TSV files
    └── mock_responses/ # Mock API responses
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

### With Markers
```bash
# Only unit tests
pytest -m unit

# Skip slow tests
pytest -m "not slow"

# Only fast tests
pytest -m "unit and not slow"
```

### With Coverage
```bash
pytest --cov=amprenta_rag --cov-report=html
open htmlcov/index.html
```

## Writing Tests

### Unit Test Example

```python
# tests/unit/ingestion/test_metabolomics_ingestion.py
import pytest
from amprenta_rag.ingestion.metabolomics.ingestion import ingest_metabolomics_file

@pytest.mark.unit
def test_ingest_metabolomics_file_success(sample_metabolomics_data):
    """Test successful metabolomics ingestion."""
    result = ingest_metabolomics_file(
        file_path=str(sample_metabolomics_data),
        create_page=False,
    )
    assert result is not None
```

### Integration Test Example

```python
# tests/integration/postgres/test_dataset_creation.py
import pytest
from amprenta_rag.ingestion.postgres_integration import create_or_update_dataset_in_postgres

@pytest.mark.integration
@pytest.mark.requires_postgres
def test_create_dataset_in_postgres(db_session):
    """Test creating dataset in Postgres."""
    dataset = create_or_update_dataset_in_postgres(
        name="Test Dataset",
        omics_type="METABOLOMICS",
        file_paths=["test.csv"],
    )
    assert dataset.id is not None
```

## Test Configuration

- `pytest.ini`: Pytest configuration
- `conftest.py`: Shared fixtures and setup
- `.env.test`: Test environment variables (not committed)

## Test Data

Sample data files are in `tests/fixtures/sample_data/`:
- `metabolomics_sample.csv`
- `proteomics_sample.csv`
- `transcriptomics_sample.csv`
- `lipidomics_sample.csv`

## Best Practices

1. **Isolation**: Each test should be independent
2. **Naming**: Use descriptive test names
3. **Fixtures**: Use fixtures for setup/teardown
4. **Mocks**: Mock external services in unit tests
5. **Markers**: Use markers to organize tests
6. **Documentation**: Add docstrings to test functions

## CI/CD

Tests run automatically on:
- Pull requests
- Commits to main branch
- Manual triggers

See `.github/workflows/test.yml` for CI configuration.

## Resources

- [Testing Guide](../docs/TESTING_GUIDE.md): Comprehensive testing documentation
- [pytest Documentation](https://docs.pytest.org/): Pytest official docs
- [Coverage.py](https://coverage.readthedocs.io/): Code coverage tool

