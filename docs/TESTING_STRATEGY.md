# Testing Strategy

Comprehensive testing approach for the Amprenta RAG System to ensure production readiness and reliability.

## Table of Contents

1. [Testing Philosophy](#testing-philosophy)
2. [Test Types](#test-types)
3. [Test Organization](#test-organization)
4. [Running Tests](#running-tests)
5. [Continuous Integration](#continuous-integration)

## Testing Philosophy

Our testing approach follows these principles:

1. **Comprehensive Coverage**: Test all major features and workflows
2. **Real-World Scenarios**: Use realistic data and conditions
3. **Production Readiness**: Ensure robustness before deployment
4. **Performance Awareness**: Measure and optimize performance
5. **Error Resilience**: Test graceful error handling

## Test Types

### Unit Tests

Test individual functions and modules in isolation.

**Location**: `amprenta_rag/tests/`

**Coverage**:
- Feature normalization functions
- Signature loading and parsing
- Metadata helpers
- Cache operations
- Utility functions

**Example**:
```bash
pytest amprenta_rag/tests/ingestion/test_feature_normalization.py -v
```

### Integration Tests

Test interactions between multiple modules and external services.

**Location**: `amprenta_rag/tests/integration/` (to be created)

**Coverage**:
- End-to-end ingestion workflows
- Notion ↔ Pinecone consistency
- Multi-module interactions
- External API integrations

**Example**:
```bash
pytest amprenta_rag/tests/integration/ -v
```

### Edge Case Tests

Test boundary conditions and unusual inputs.

**Location**: `amprenta_rag/tests/edge_cases/` (to be created)

**Coverage**:
- Empty files
- Malformed data
- Extreme values
- Missing optional fields
- Invalid API responses

**Example**:
```bash
pytest amprenta_rag/tests/edge_cases/ -v
```

### Performance Tests

Measure and benchmark performance characteristics.

**Location**: `amprenta_rag/tests/performance/` (to be created)

**Coverage**:
- Feature cache performance
- Batch ingestion throughput
- Query response times
- Memory footprint
- Concurrent operation handling

**Example**:
```bash
pytest amprenta_rag/tests/performance/ -v --benchmark-only
```

### Load/Stress Tests

Test system behavior under high load.

**Location**: `scripts/load_tests/` (to be created)

**Coverage**:
- Large batch operations
- Concurrent requests
- Memory usage under load
- API rate limit handling
- Recovery from failures

**Example**:
```bash
python scripts/load_tests/test_batch_ingestion_load.py
```

## Test Organization

```
amprenta_rag/tests/
├── __init__.py
├── ingestion/
│   ├── test_feature_extraction.py
│   ├── test_feature_normalization.py
│   ├── test_dataset_feature_cache.py
│   └── ...
├── integration/          # NEW
│   ├── test_end_to_end_ingestion.py
│   ├── test_notion_pinecone_sync.py
│   └── ...
├── edge_cases/          # NEW
│   ├── test_empty_files.py
│   ├── test_malformed_data.py
│   └── ...
├── performance/         # NEW
│   ├── test_cache_performance.py
│   ├── test_batch_throughput.py
│   └── ...
└── fixtures/            # NEW
    ├── sample_data.py
    └── mock_services.py
```

## Running Tests

### Run All Tests

```bash
# Run all tests
pytest amprenta_rag/tests/ -v

# Run with coverage
pytest --cov=amprenta_rag --cov-report=html

# Run specific test file
pytest amprenta_rag/tests/ingestion/test_feature_extraction.py -v

# Run specific test function
pytest amprenta_rag/tests/ingestion/test_feature_extraction.py::test_extract_lipidomics -v
```

### Run Integration Tests

```bash
# All integration tests
pytest amprenta_rag/tests/integration/ -v

# Integration tests with real services (requires API keys)
pytest amprenta_rag/tests/integration/ --real-services -v
```

### Run Performance Tests

```bash
# Performance benchmarks
pytest amprenta_rag/tests/performance/ -v --benchmark-only

# Compare performance
pytest amprenta_rag/tests/performance/ -v --benchmark-compare
```

### Run Comprehensive Feature Test

```bash
# Test all features with real data
python scripts/test_all_features.py
```

## Test Categories

### Critical Path Tests

These tests must pass before any deployment:

1. ✅ Feature extraction from datasets
2. ✅ Signature scoring accuracy
3. ✅ Multi-omics normalization
4. ✅ Cache operations
5. ✅ Basic RAG queries

### Extended Tests

Run periodically or before major releases:

1. Load/stress testing
2. Performance benchmarks
3. Edge case coverage
4. Integration workflows
5. Error recovery scenarios

### Development Tests

Run during development:

1. Unit tests for new features
2. Quick integration checks
3. Code quality checks

## Continuous Integration

### Test Pipeline

```yaml
# .github/workflows/tests.yml (to be created)
- Unit tests (fast, always run)
- Integration tests (requires API keys, run on PR)
- Performance benchmarks (run nightly)
- Load tests (run weekly)
```

### Pre-commit Checks

```bash
# Run quick checks before commit
pytest amprenta_rag/tests/ -k "not integration" --fast
```

## Test Data Management

### Synthetic Test Data

Use synthetic data for unit tests to avoid external dependencies:

```python
# fixtures/sample_data.py
SAMPLE_LIPIDOMICS_CSV = """
Species,Intensity
Cer(d18:1/16:0),1000
Cer(d18:1/18:0),1500
"""

SAMPLE_SIGNATURE_TSV = """
feature_type	feature_name	direction	weight
lipid	Cer(d18:1/16:0)	↑	1.0
"""
```

### Real Test Data

Use real data for integration tests (with proper cleanup):

- Mark test data with special tags
- Clean up after tests
- Use dedicated test databases when possible

## Performance Benchmarks

### Key Metrics

1. **Feature Extraction Speed**
   - Single dataset: < 5 seconds
   - Batch (10 datasets): < 30 seconds

2. **Signature Scoring Speed**
   - Single signature: < 2 seconds
   - Batch (100 signatures): < 60 seconds

3. **Query Response Time**
   - Simple query: < 3 seconds
   - Complex query with filters: < 10 seconds

4. **Memory Usage**
   - Feature cache: < 500MB for 1000 datasets
   - Batch ingestion: < 2GB for 100 files

### Benchmark Targets

| Operation | Target | Current |
|-----------|--------|---------|
| Feature extraction | < 5s/dataset | TBD |
| Signature scoring | < 2s/signature | TBD |
| Batch ingestion | < 1min/10 files | TBD |
| Query response | < 3s | TBD |

## Error Scenario Testing

### Network Failures

- Simulate API timeouts
- Test retry logic
- Verify graceful degradation

### Invalid Data

- Malformed CSV/TSV files
- Missing required fields
- Invalid feature names
- Corrupted cache files

### Configuration Errors

- Missing API keys
- Invalid database IDs
- Malformed configuration

## Test Coverage Goals

| Category | Current | Target |
|----------|---------|--------|
| Overall Coverage | 39% | 80%+ |
| Unit Tests | ~40% | 80%+ |
| Integration Tests | ~20% | 60%+ |
| Edge Cases | ~10% | 50%+ |
| Performance Tests | 0% | 100% |
| New Modules | - | 80%+ (required) |

**CI Enforcement**: 35% threshold (rising incrementally toward 80%)

## Best Practices

1. **Write tests first** for new features (TDD required for new features)
2. **Keep tests fast** - unit tests should run in < 1 second each
3. **Use fixtures** for common test data
4. **Mock external services** in unit tests
5. **Clean up** test data after tests
6. **Document** test purpose and assumptions
7. **Run tests frequently** during development
8. **Minimum 80% coverage for new modules**
9. **CI enforces coverage threshold** (35% baseline, rising incrementally)

## Next Steps

1. Create test infrastructure for new test types
2. Expand edge case coverage
3. Implement performance benchmarks
4. Add load/stress testing
5. Set up CI/CD pipeline

## Resources

- [pytest Documentation](https://docs.pytest.org/)
- [Testing Best Practices](https://docs.python-guide.org/writing/tests/)
- [Performance Testing Guide](docs/PERFORMANCE_TESTING.md) (to be created)

