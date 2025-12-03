# Refactoring, Cleanup, Testing & Documentation Plan

**Date**: December 2, 2025  
**Current State**: Feature-complete, production-ready, but needs polish  
**Goal**: Improve code quality, test coverage, documentation, and maintainability

---

## 📊 **CURRENT STATE ASSESSMENT**

### Codebase Statistics:
- **Python Files**: 38 files
- **Total Lines**: ~9,400 lines
- **Test Files**: 3 tests (minimal coverage)
- **Functions/Classes**: 93+ in ingestion module alone
- **Documentation**: Inconsistent (some modules well-documented, others sparse)

### Test Coverage:
- ✅ **3 tests passing**: `test_sanitize_metadata`, `test_build_meta_filter` (2 tests)
- ❌ **Coverage**: Very low (~5-10% estimated)
- ❌ **Missing**: Tests for ingestion pipelines, signature matching, scoring, etc.

### Documentation Status:
- ✅ **Some modules**: Good docstrings (e.g., `signature_matching.py`)
- ⚠️ **Many modules**: Minimal or missing docstrings
- ⚠️ **README**: Exists but could be more comprehensive
- ❌ **API Documentation**: Missing

### Code Quality:
- ✅ **Structure**: Well-organized into modules
- ⚠️ **Duplication**: Some shared logic could be extracted
- ⚠️ **Error Handling**: Inconsistent patterns
- ⚠️ **Type Hints**: Partial (some functions have them, others don't)

---

## 🎯 **REFACTORING PRIORITIES**

### Phase 1: Code Cleanup & Organization (High Priority)

#### 1.1 Remove Dead Code
- [ ] Search for unused imports
- [ ] Remove commented-out code blocks
- [ ] Remove unused functions/classes
- [ ] Clean up temporary test files

#### 1.2 Standardize Code Style
- [ ] Run black formatter across all files
- [ ] Fix line length issues
- [ ] Standardize import ordering (isort)
- [ ] Consistent naming conventions

#### 1.3 Extract Common Patterns
- [ ] **Notion API calls**: Create shared wrapper for common patterns
- [ ] **Error handling**: Standardize try/except patterns
- [ ] **Logging**: Ensure consistent logging patterns
- [ ] **Batch processing**: Extract batching logic to shared utility

#### 1.4 Type Hints
- [ ] Add type hints to all public functions
- [ ] Add type hints to module-level functions
- [ ] Use `typing` module consistently
- [ ] Add return type annotations

---

### Phase 2: Documentation (High Priority)

#### 2.1 Module-Level Docstrings
- [ ] Add docstrings to all modules explaining purpose
- [ ] Document module-level constants
- [ ] Add usage examples where appropriate

#### 2.2 Function/Class Docstrings
- [ ] Add docstrings to all public functions
- [ ] Add docstrings to all classes
- [ ] Use Google/NumPy docstring style consistently
- [ ] Document parameters, return values, exceptions

#### 2.3 Inline Comments
- [ ] Add comments to complex logic
- [ ] Explain "why" not just "what"
- [ ] Document algorithm choices
- [ ] Add TODO comments for future improvements

#### 2.4 API Documentation
- [ ] Create API reference documentation
- [ ] Document public interfaces
- [ ] Add usage examples
- [ ] Document configuration options

---

### Phase 3: Testing (High Priority)

#### 3.1 Unit Tests
- [ ] **Ingestion Modules**:
  - [ ] `dataset_ingestion.py`: Test mwTab extraction, chunking, embedding
  - [ ] `zotero_ingest.py`: Test attachment processing, metadata extraction
  - [ ] `email_ingestion.py`: Test email parsing, chunking
  - [ ] `experiments_ingestion.py`: Test experiment data extraction
  - [ ] `signature_ingestion.py`: Test signature loading, component creation
  - [ ] `signature_matching.py`: Test matching logic, scoring
  - [ ] `signature_detection.py`: Test detection algorithms
  - [ ] `feature_extraction.py`: Test metabolite extraction

- [ ] **Query Modules**:
  - [ ] `rag_engine.py`: Test query orchestration
  - [ ] `pinecone_query.py`: Test filtering, embedding

- [ ] **Signature Modules**:
  - [ ] `signature_scoring.py`: Test scoring algorithms
  - [ ] `species_matching.py`: Test matching logic
  - [ ] `signature_loader.py`: Test loading from various formats

- [ ] **Utility Modules**:
  - [ ] `metadata_semantic.py`: Test metadata extraction
  - [ ] `notion_pages.py`: Test Notion API helpers
  - [ ] `pinecone_utils.py`: Test metadata sanitization
  - [ ] `text_extraction.py`: Test text extraction

#### 3.2 Integration Tests
- [ ] End-to-end ingestion pipeline tests
- [ ] Signature matching integration tests
- [ ] RAG query integration tests
- [ ] Notion API integration tests (mocked)

#### 3.3 Test Infrastructure
- [ ] Set up pytest fixtures for common test data
- [ ] Create mock Notion API responses
- [ ] Create mock Pinecone responses
- [ ] Set up test data files
- [ ] Configure test coverage reporting

---

### Phase 4: Code Refactoring (Medium Priority)

#### 4.1 Extract Shared Utilities
- [ ] **Notion API Wrapper**: Centralize all Notion API calls
- [ ] **Pinecone Operations**: Centralize Pinecone operations
- [ ] **Batch Processing**: Extract batching logic
- [ ] **Error Handling**: Standardize error handling patterns

#### 4.2 Reduce Duplication
- [ ] **Text Chunking**: Already shared, verify consistency
- [ ] **Embedding**: Already shared, verify consistency
- [ ] **Metadata Extraction**: Check for duplication
- [ ] **Notion Page Creation**: Standardize patterns

#### 4.3 Improve Error Handling
- [ ] Consistent exception types
- [ ] Better error messages
- [ ] Proper exception chaining
- [ ] Graceful degradation where appropriate

#### 4.4 Configuration Management
- [ ] Validate configuration on startup
- [ ] Better error messages for missing config
- [ ] Configuration documentation
- [ ] Environment variable validation

---

### Phase 5: Performance & Optimization (Low Priority)

#### 5.1 Performance Improvements
- [ ] Profile slow operations
- [ ] Optimize batch sizes
- [ ] Cache expensive operations
- [ ] Parallel processing where appropriate

#### 5.2 Resource Management
- [ ] Connection pooling for APIs
- [ ] Rate limiting handling
- [ ] Memory optimization for large datasets

---

## 📋 **DETAILED TASK BREAKDOWN**

### Immediate Actions (This Week)

#### Day 1: Code Cleanup
1. Run linters (flake8, pylint, mypy)
2. Fix all linting errors
3. Remove unused imports
4. Remove commented code
5. Run black formatter

#### Day 2: Documentation
1. Add module docstrings to all ingestion modules
2. Add function docstrings to public APIs
3. Document configuration options
4. Update README with examples

#### Day 3: Testing Foundation
1. Set up pytest configuration
2. Create test fixtures
3. Create mock helpers
4. Write tests for core utilities (5-10 tests)

#### Day 4-5: Core Tests
1. Write tests for signature matching (5-10 tests)
2. Write tests for metadata extraction (5-10 tests)
3. Write tests for text processing (3-5 tests)
4. Write integration tests (2-3 tests)

### Short-term (Next 2 Weeks)

1. Complete unit test coverage for all ingestion modules
2. Add integration tests for end-to-end workflows
3. Complete API documentation
4. Refactor shared utilities
5. Improve error handling

### Long-term (Next Month)

1. Achieve 70%+ test coverage
2. Complete all documentation
3. Performance optimization
4. Advanced refactoring

---

## 🛠️ **TOOLS & SETUP**

### Linting & Formatting:
```bash
# Install tools
pip install black isort flake8 pylint mypy pytest pytest-cov

# Format code
black amprenta_rag/ scripts/
isort amprenta_rag/ scripts/

# Lint
flake8 amprenta_rag/ scripts/
pylint amprenta_rag/

# Type checking
mypy amprenta_rag/
```

### Testing:
```bash
# Run tests
pytest amprenta_rag/tests -v

# With coverage
pytest amprenta_rag/tests --cov=amprenta_rag --cov-report=html

# Specific test
pytest amprenta_rag/tests/test_signature_matching.py -v
```

### Documentation:
```bash
# Generate API docs (if using Sphinx)
sphinx-build -b html docs/ docs/_build/
```

---

## 📊 **SUCCESS METRICS**

### Code Quality:
- ✅ Zero linting errors
- ✅ 100% type hint coverage for public APIs
- ✅ Consistent code style (black-formatted)
- ✅ No duplicate code patterns

### Documentation:
- ✅ All modules have docstrings
- ✅ All public functions have docstrings
- ✅ README is comprehensive
- ✅ API documentation exists

### Testing:
- ✅ 70%+ test coverage
- ✅ All critical paths tested
- ✅ Integration tests for main workflows
- ✅ Tests run in < 30 seconds

### Maintainability:
- ✅ Clear module boundaries
- ✅ Minimal code duplication
- ✅ Consistent error handling
- ✅ Easy to add new features

---

## 🎯 **RECOMMENDED STARTING POINT**

### Option A: Quick Wins (Recommended)
Start with **Phase 1.1-1.2** (Code Cleanup):
1. Run black formatter
2. Remove unused imports
3. Fix linting errors
4. Add basic docstrings

**Time**: 2-4 hours  
**Impact**: Immediate code quality improvement

### Option B: Testing First
Start with **Phase 3.1** (Unit Tests):
1. Set up test infrastructure
2. Write tests for core utilities
3. Write tests for signature matching
4. Achieve 30%+ coverage

**Time**: 1-2 days  
**Impact**: Prevents regressions, enables confident refactoring

### Option C: Documentation First
Start with **Phase 2** (Documentation):
1. Add module docstrings
2. Add function docstrings
3. Update README
4. Document configuration

**Time**: 1-2 days  
**Impact**: Improves developer experience, onboarding

---

## 💡 **RECOMMENDATION**

**Start with Option A (Quick Wins)**:
1. Run formatters and fix style issues (1 hour)
2. Add basic docstrings to key modules (2-3 hours)
3. Then move to testing (Option B)

This gives immediate visible improvements while setting up for more substantial work.

---

**Status**: Ready to begin  
**Estimated Total Time**: 1-2 weeks for comprehensive cleanup  
**Priority**: High (improves maintainability and developer experience)

