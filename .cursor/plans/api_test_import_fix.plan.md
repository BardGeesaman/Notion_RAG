# API Test Import Chain Fix

## Problem
When `TestClient(app)` imports `amprenta_rag.api.main`, it triggers import of ALL 90+ routers. Some routers have module-level imports that cascade into:
- Celery tasks (require broker connection)
- Database sessions (attempt connection)
- ML model loading (heavy dependencies)
- External service clients

This breaks API tests even though the tests only need 1-2 routers.

**Affected tests:**
- `test_retrosynthesis_api.py`
- `test_sar_api.py`
- `test_gnn_api.py`
- `test_active_learning_api.py`

---

## Solution: App Factory Pattern

### 1. Create App Factory
Convert `main.py` to use factory pattern:

```python
# amprenta_rag/api/main.py

from typing import Optional, List, Set

# Router registry (lazy import)
ROUTER_REGISTRY = {
    "admin": ("amprenta_rag.api.routers.admin", "router", "/api/v1", ["Admin"]),
    "programs": ("amprenta_rag.api.routers.programs", "router", "/api/v1/programs", ["Programs"]),
    "retrosynthesis": ("amprenta_rag.api.routers.retrosynthesis", "router", "/api/v1/retrosynthesis", ["Retrosynthesis"]),
    # ... all 90+ routers
}

def create_app(
    include_routers: Optional[Set[str]] = None,
    exclude_routers: Optional[Set[str]] = None,
    skip_middleware: bool = False,
) -> FastAPI:
    """Create FastAPI app with configurable routers.
    
    Args:
        include_routers: If provided, ONLY include these routers
        exclude_routers: Exclude these routers (ignored if include_routers set)
        skip_middleware: Skip rate limiting and other middleware (for tests)
    
    Returns:
        Configured FastAPI application
    """
    app = FastAPI(
        title="Amprenta Multi-Omics Platform API",
        version="1.0.0",
    )
    
    # Add middleware (optional)
    if not skip_middleware:
        app.add_middleware(RequestSizeLimitMiddleware)
        app.add_middleware(TimeoutMiddleware)
        # ... etc
    
    # Determine which routers to include
    if include_routers:
        router_names = include_routers
    else:
        router_names = set(ROUTER_REGISTRY.keys())
        if exclude_routers:
            router_names -= exclude_routers
    
    # Lazy import and register routers
    for name in router_names:
        if name not in ROUTER_REGISTRY:
            continue
        
        module_path, router_attr, prefix, tags = ROUTER_REGISTRY[name]
        try:
            module = importlib.import_module(module_path)
            router = getattr(module, router_attr)
            app.include_router(router, prefix=prefix, tags=tags)
        except ImportError as e:
            logger.warning(f"Failed to import router {name}: {e}")
    
    return app

# Default app for production (includes all routers)
app = create_app()
```

### 2. Create Test Helper
```python
# amprenta_rag/tests/conftest.py

import pytest
from fastapi.testclient import TestClient

def create_test_app(routers: List[str]) -> FastAPI:
    """Create minimal test app with only specified routers."""
    from amprenta_rag.api.main import create_app
    return create_app(
        include_routers=set(routers),
        skip_middleware=True
    )

@pytest.fixture
def retrosynthesis_client():
    """Client with only retrosynthesis router."""
    app = create_test_app(["retrosynthesis"])
    # Mock auth
    app.dependency_overrides[get_current_user] = mock_user
    return TestClient(app)

@pytest.fixture  
def sar_client():
    """Client with only SAR router."""
    app = create_test_app(["sar"])
    app.dependency_overrides[get_current_user] = mock_user
    return TestClient(app)
```

### 3. Update Affected Tests
```python
# test_retrosynthesis_api.py

def test_analyze_target(retrosynthesis_client):
    """Test with minimal app - no import chain issues."""
    response = retrosynthesis_client.post(
        "/api/v1/retrosynthesis/analyze",
        json={"smiles": "CCO"}
    )
    assert response.status_code == 200
```

---

## Implementation Batches

### Batch 1: Create Router Registry (Day 1)
- [ ] Extract all router registrations to `ROUTER_REGISTRY` dict
- [ ] Create `create_app()` function with include/exclude support
- [ ] Verify production app still works (`uvicorn amprenta_rag.api.main:app`)

### Batch 2: Test Infrastructure (Day 2)
- [ ] Create `create_test_app()` helper in conftest.py
- [ ] Create router-specific fixtures (retrosynthesis_client, sar_client, etc.)
- [ ] Add mock user fixture

### Batch 3: Fix Affected Tests (Day 3)
- [ ] Update test_retrosynthesis_api.py to use minimal app
- [ ] Update test_sar_api.py
- [ ] Update test_gnn_api.py
- [ ] Update test_active_learning_api.py
- [ ] Run all API tests

### Batch 4: Verify & Document (Day 4)
- [ ] Run full test suite
- [ ] Update testing documentation
- [ ] Remove tech debt item from ROADMAP

---

## Success Criteria

- [ ] All previously blocked API tests pass
- [ ] Production app unchanged (all routers load)
- [ ] Test isolation (each test only loads required routers)
- [ ] Test speed improvement (fewer imports per test)

---

## Alternative Approaches (Considered)

### A. Lazy Imports in Every Router
**Pros:** No app changes
**Cons:** 90+ routers to modify, easy to miss one, ongoing maintenance

### B. TESTING Environment Variable
**Pros:** Simple
**Cons:** Still imports all modules, just skips some initialization

### C. Separate Test App Module
**Pros:** Clean separation
**Cons:** Duplicate registration code, drift risk

**Chosen: App Factory (Option above)** - Cleanest, most flexible, standard pattern.

---

## TODOs

- [ ] Batch 1: Router Registry + create_app()
- [ ] Batch 2: Test Infrastructure  
- [ ] Batch 3: Fix Affected Tests
- [ ] Batch 4: Verify & Document


