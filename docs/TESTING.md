# Testing Documentation

This document captures testing patterns, best practices, and future improvements for the Amprenta RAG system.

---

## E2E Test Pattern: In-Memory API Stub

For Streamlit E2E tests needing API data without DB dependencies:

1. Create session-scoped FastAPI fixture with dict storage
2. Define minimal CRUD endpoints matching API schema
3. Run uvicorn in daemon thread on separate port
4. Set Streamlit's `API_URL` env var to stub URL
5. Use `_kill_port()` cleanup before starting
6. Health check polling before yielding

**Reference:** `amprenta_rag/tests/dashboard/test_e2e_dvc.py`

**Benefits:** Zero dependencies, fast (<1s), CI-friendly, isolated

### Future Improvements (P3)

1. **Schema Validation Test** - Compare stub response to DB model schema
2. **Extract to Shared Fixture** - `tests/fixtures/api_stub.py` for reuse
3. **Add Endpoints as Needed** - Keep stub minimal, extend when tests require

---

## DVC Integration - Future Work (P3)

1. **Add Unit Tests** - Mock subprocess for dvc_manager functions
2. **Add Index** - `CREATE INDEX ix_datasets_dvc_version ON datasets(dvc_version)`
3. **Progress Spinner** - Show spinner during DVC restore
4. **Dataset Diff Viewer** - Compare two dataset versions (Phase 3)
5. **S3 Remote Storage** - Phase 2 enhancement

---

For comprehensive testing strategy and guidelines, see:
- [TESTING_STRATEGY.md](TESTING_STRATEGY.md) - Overall testing approach and TDD requirements
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Detailed testing instructions and examples

