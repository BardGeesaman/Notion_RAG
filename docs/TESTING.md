# Testing Best Practices

## Table of Contents
- [Environment Setup](#environment-setup)
- [Test Organization](#test-organization)
- [Database Testing](#database-testing)
- [API Testing](#api-testing)
- [E2E Testing](#e2e-testing)
- [Common Issues & Fixes](#common-issues--fixes)

---

## Environment Setup

### Automatic Environment Activation with direnv

This project uses `direnv` to automatically activate the conda environment when you enter the project directory.

**One-time setup**:

1. Install direnv:
   ```bash
   # macOS
   brew install direnv
   
   # Ubuntu/Debian
   sudo apt-get install direnv
   
   # Other: see https://direnv.net/docs/installation.html
   ```

2. Add direnv hook to your shell RC file:
   ```bash
   # For zsh (add to ~/.zshrc)
   eval "$(direnv hook zsh)"
   
   # For bash (add to ~/.bashrc)
   eval "$(direnv hook bash)"
   
   # For fish (add to ~/.config/fish/config.fish)
   direnv hook fish | source
   ```

3. Restart your shell or source the RC file:
   ```bash
   source ~/.zshrc  # or ~/.bashrc
   ```

4. Allow direnv in the project directory:
   ```bash
   cd /path/to/RAG
   direnv allow
   ```

**Verification**:
```bash
# cd into project - should auto-activate
cd /path/to/RAG
# You should see: direnv: loading .envrc
# Prompt should show: (myenv)

# cd out of project - should auto-deactivate
cd ~
# You should see: direnv: unloading
```

**Troubleshooting**:

- **"direnv: error .envrc is blocked"**: Run `direnv allow` in the project directory
- **"conda: command not found"**: Ensure conda is installed and initialized in your shell RC file
- **Environment doesn't activate**: Check that direnv hook is in your shell RC file and shell is restarted
- **"myenv not found"**: Create the conda environment first (see below)

**Creating the conda environment**:

```bash
# Using environment.yml (recommended)
conda env create -f environment.yml

# OR manually create environment
conda create -n myenv python=3.11
conda activate myenv
pip install -r requirements.txt
```

### Manual Activation (if not using direnv)

```bash
# Activate conda environment
conda activate myenv

# Install test dependencies
pip install pytest pytest-asyncio playwright httpx

# For E2E tests with Playwright
playwright install
```

### Run Tests
```bash
# Run all tests
pytest

# Run specific test file
pytest amprenta_rag/tests/services/test_activity_service.py

# Run with verbose output
pytest -v

# Run with print statements visible
pytest -s

# Run single test
pytest path/to/test.py::TestClass::test_method
```

---

## Test Organization

### Module Naming
**CRITICAL**: Test module names must be unique across the entire test suite.

❌ **BAD**:
```
amprenta_rag/tests/utils/test_activity.py
amprenta_rag/tests/services/test_activity.py  # Duplicate basename!
```

✅ **GOOD**:
```
amprenta_rag/tests/utils/test_activity_utils.py
amprenta_rag/tests/services/test_activity_service.py
```

**Why**: Python's import system gets confused by duplicate basenames, causing `ImportError` with `__file__ mismatch`.

**Fix**: Rename conflicting files and clear cache:
```bash
find amprenta_rag/tests -name "__pycache__" -exec rm -rf {} +
```

---

## Database Testing

### Foreign Key Constraints
Foreign keys require **real records** or `None`, not random UUIDs.

❌ **BAD**:
```python
log_activity(
    actor_id=uuid4(),  # FK violation if user doesn't exist!
    program_id=uuid4(),  # FK violation if program doesn't exist!
)
```

✅ **GOOD**:
```python
# Option 1: Use None for system actions
log_activity(
    actor_id=None,
    program_id=None,
)

# Option 2: Create real records first
user = User(id=uuid4(), username="test_user", ...)
db.add(user)
db.commit()

log_activity(actor_id=user.id, ...)
```

### Test Data Uniqueness
Always use UUID-based unique values for constrained fields.

❌ **BAD**:
```python
user = User(username="test_user", email="test@test.com")  # Fails on 2nd run!
```

✅ **GOOD**:
```python
from uuid import uuid4

user = User(
    username=f"testuser_{uuid4().hex[:8]}",
    email=f"test_{uuid4().hex[:8]}@test.com",
)
```

### Session Management
Objects accessed outside their session must be detached with `db.expunge()`.

❌ **BAD**:
```python
def get_items():
    with db_session() as db:
        items = db.query(Item).all()
        return items  # DetachedInstanceError when accessed!
```

✅ **GOOD**:
```python
def get_items():
    with db_session() as db:
        items = db.query(Item).all()
        for item in items:
            db.expunge(item)  # Detach before session closes
        return items
```

### SQLAlchemy Comparisons
Use `.is_(False)` for boolean comparisons, not `== False`.

❌ **BAD**:
```python
query.filter(Model.is_active == False)  # Linting error E712
```

✅ **GOOD**:
```python
query.filter(Model.is_active.is_(False))
```

### Complex Updates with Joins
`.update()` doesn't work reliably with joins. Use two-step approach:

❌ **BAD**:
```python
count = (
    db.query(Model1)
    .join(Model2)
    .filter(Model2.user_id == user_id)
    .update({"is_read": True})  # Returns 0!
)
```

✅ **GOOD**:
```python
# Step 1: Fetch IDs
ids = (
    db.query(Model1.id)
    .join(Model2)
    .filter(Model2.user_id == user_id)
    .all()
)

# Step 2: Update by ID
ids_list = [id[0] for id in ids]
count = (
    db.query(Model1)
    .filter(Model1.id.in_(ids_list))
    .update({"is_read": True}, synchronize_session=False)
)
```

---

## API Testing

### Authentication with Dependency Overrides
Override `get_current_user` dependency for tests.

```python
from fastapi.testclient import TestClient
from amprenta_rag.api.main import app
from amprenta_rag.api.dependencies import get_current_user

# Mock user
TEST_USER_ID = uuid4()

def mock_current_user():
    class FakeUser:
        id = TEST_USER_ID
        username = "test_user"
        email = "test@example.com"
    return FakeUser()

def test_authenticated_endpoint():
    # Override dependency
    app.dependency_overrides[get_current_user] = mock_current_user
    
    try:
        client = TestClient(app)
        response = client.get("/api/v1/protected")
        assert response.status_code == 200
    finally:
        # Clean up
        app.dependency_overrides.clear()
```

### Database Mocking
Mock `get_db` dependency to avoid real DB hits in unit tests.

```python
from amprenta_rag.database.base import get_db

def test_with_mock_db():
    mock_db = MagicMock()
    
    # Configure mock query chain
    mock_db.query.return_value.filter.return_value.first.return_value = mock_user
    
    # Override dependency
    app.dependency_overrides[get_db] = lambda: mock_db
    
    try:
        client = TestClient(app)
        response = client.post("/api/v1/endpoint")
        assert response.status_code == 200
    finally:
        app.dependency_overrides.clear()
```

---

## E2E Testing

### Selector Best Practices
Use semantic selectors, not implementation details.

❌ **BAD**:
```python
page.locator('input[aria-label="Search"]')  # Fragile, implementation-specific
page.locator('#submit-btn')  # ID can change
```

✅ **GOOD**:
```python
page.get_by_text("Search")  # Semantic, user-facing
page.get_by_role("button", name="Submit")  # ARIA role
page.get_by_placeholder("Enter email")  # Placeholder text
```

### Handling Optional Elements
Use `.or_()` pattern for elements that may or may not exist.

```python
# Check for either success message or error
page.locator('text="Success"').or_(page.locator('text="Error"')).wait_for()

# Or use try/except
try:
    page.get_by_text("Expected element", timeout=5000).wait_for()
except TimeoutError:
    # Element doesn't exist, handle gracefully
    pass
```

### Scoping Selectors
Scope to main content to avoid matching hidden sidebar elements.

❌ **BAD**:
```python
page.get_by_role("tab", name="Interactive")  # Might match hidden sidebar!
```

✅ **GOOD**:
```python
main = page.locator('[data-testid="stMainBlockContainer"]').first
tab = main.get_by_role("tab", name="Interactive")
```

### Wait Strategies
Use appropriate waits for different scenarios.

```python
# Basic page load
page.goto(url)
page.wait_for_load_state("domcontentloaded")

# Wait for specific element
page.wait_for_selector('[data-testid="stAppViewContainer"]', timeout=15000)

# networkidle is often too slow/fragile - avoid unless necessary
# page.wait_for_load_state("networkidle")  # Avoid!
```

### Streamlit-Specific Patterns

#### Navigation with Query Parameters
```python
# Direct navigation to specific page
page.goto(f"{base_url}/?page=Datasets")
page.wait_for_load_state("domcontentloaded")
```

#### Triggering Reruns
```python
# After filling textarea, press Tab to trigger Streamlit rerun
page.fill('textarea[aria-label="Comment"]', "Test comment")
page.keyboard.press("Tab")  # Triggers Streamlit to process input
```

#### Helper Pattern
```python
def _main_container(page):
    """Scope selectors to main content area."""
    return page.locator('[data-testid="stMainBlockContainer"]').first

# Usage
main = _main_container(page)
button = main.get_by_role("button", name="Submit")
```

---

## Common Issues & Fixes

### Import Collisions
**Error**: `ImportError: imported module 'test_foo' has __file__ mismatch`

**Cause**: Duplicate test module basenames in different directories.

**Fix**:
```bash
# 1. Clear all __pycache__
find amprenta_rag/tests -name "__pycache__" -exec rm -rf {} +

# 2. Rename conflicting files
mv test_activity.py test_activity_service.py

# 3. Run tests again
pytest
```

### DetachedInstanceError
**Error**: `Instance <Model at 0x...> is not bound to a Session`

**Cause**: Accessing SQLAlchemy object outside its session context.

**Fix**: Use `db.expunge(obj)` before returning from session.

### Foreign Key Violations
**Error**: `ForeignKeyViolation: Key (field_id)=(uuid) is not present in table "table"`

**Cause**: Using `uuid4()` for FK field without creating the referenced record.

**Fix**: Either use `None` or create the referenced record first.

### Unique Constraint Violations
**Error**: `IntegrityError: duplicate key value violates unique constraint`

**Cause**: Hardcoded test data (e.g., `username="test_user"`) used in multiple tests.

**Fix**: Use UUID-based unique values:
```python
username = f"testuser_{uuid4().hex[:8]}"
```

### Linting: Unused Variables
**Error**: `F841 Local variable 'foo' is assigned to but never used`

**Fix**: Don't assign if you don't use the result:
```python
# Before
result = some_function()  # Unused

# After
some_function()
```

### Boolean Comparisons
**Error**: `E712 Avoid equality comparisons to False`

**Fix**: Use `.is_(False)`:
```python
# Before
filter(Model.active == False)

# After
filter(Model.active.is_(False))
```

### Mock Import Errors
**Error**: `AttributeError: module 'pytest' has no attribute 'mock'`

**Cause**: Using `pytest.mock` which doesn't exist.

**Fix**: Use unittest.mock:
```python
# WRONG
import pytest.mock
pytest.mock.patch(...)

# CORRECT
from unittest.mock import patch, MagicMock, ANY

@patch("module.function")
def test_foo(mock_func):
    mock_func.return_value = "test"
```

### Mock Patch Path Errors
**Error**: Mock doesn't intercept calls, assertions fail

**Cause**: Patching where function is DEFINED instead of where it's USED.

**Fix**: Patch at the import location:
```python
# Function defined in: amprenta_rag.utils.helpers
# Function imported in: amprenta_rag.services.foo

# WRONG - patches the definition
@patch("amprenta_rag.utils.helpers.my_function")

# CORRECT - patches where it's used
@patch("amprenta_rag.services.foo.my_function")
```

### SQLAlchemy Mock Chain Setup
**Error**: `db.add()` or `db.commit()` not called, mock returns MagicMock instead of data

**Cause**: SQLAlchemy query chains need full mock setup.

**Fix**: Configure the full chain:
```python
# Setup mock for: db.query(Model).filter(...).first()
mock_db.query.return_value.filter.return_value.first.return_value = mock_entity

# Setup mock for: db.query(Model).filter(...).all()
mock_db.query.return_value.filter.return_value.all.return_value = [mock_entity]

# For new record creation, return None to trigger db.add()
mock_db.query.return_value.filter.return_value.first.return_value = None
```

### Mock Return vs Side Effect
**Error**: Mock returns wrong value or doesn't raise exception

**Fix**: Use correct mock attribute:
```python
# Single return value
mock_func.return_value = "value"

# Raise exception
mock_func.side_effect = ValueError("error")

# Multiple calls with different values
mock_func.side_effect = ["first", "second", ValueError("third")]
```

---

## Test Checklist

Before committing tests:

- [ ] Test module name is unique (no duplicate basenames)
- [ ] All test data uses UUID-based unique values
- [ ] Foreign keys use `None` or real records (not random UUIDs)
- [ ] SQLAlchemy objects are detached with `db.expunge()` if accessed outside session
- [ ] API tests use `app.dependency_overrides` with `try/finally` cleanup
- [ ] E2E selectors use semantic locators (not CSS/XPath)
- [ ] E2E waits use `domcontentloaded` + specific selectors (not `networkidle`)
- [ ] No unused imports or variables
- [ ] Boolean comparisons use `.is_(False)` not `== False`
- [ ] Tests pass: `pytest path/to/test.py -v`
- [ ] Linting clean: `ruff check path/to/test.py`

---

## Running the Full Test Suite

```bash
# Full test suite
pytest amprenta_rag/tests/ -v

# Skip slow E2E tests
pytest amprenta_rag/tests/ -v -m "not e2e"

# Run only unit tests
pytest amprenta_rag/tests/unit/ -v

# Run with coverage
pytest --cov=amprenta_rag --cov-report=html

# Parallel execution (faster)
pytest -n auto
```

---

## Further Reading

- [Pytest Documentation](https://docs.pytest.org/)
- [Playwright Best Practices](https://playwright.dev/python/docs/best-practices)
- [SQLAlchemy Session Management](https://docs.sqlalchemy.org/en/20/orm/session_basics.html)
- [FastAPI Testing](https://fastapi.tiangolo.com/tutorial/testing/)
