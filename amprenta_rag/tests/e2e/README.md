# E2E Testing Patterns

This document describes patterns and conventions for Playwright E2E tests in Amprenta RAG.

## Quick Start

```bash
# Start servers
cd /Users/bard/Documents/RAG && conda activate myenv
uvicorn amprenta_rag.api.main:app --port 8000 &
DISABLE_AUTH=true streamlit run scripts/dashboard/app.py --server.port 8501 &

# Run E2E tests (headed mode for debugging)
pytest amprenta_rag/tests/e2e/ -v --headed

# Run specific test
pytest amprenta_rag/tests/e2e/test_retrosynthesis_e2e.py -v --headed
```

## Selector Patterns

### ✅ DO: Use Role-Based Selectors

```python
# Best: Role + accessible name
page.get_by_role("button", name="Submit")
page.get_by_role("heading", name="Dashboard")
page.get_by_role("tab", name="Settings")

# Good: Label-based
page.get_by_label("Email address")
page.get_by_placeholder("Enter SMILES...")

# Good: Test ID (if available)
page.get_by_test_id("submit-button")
```

### ❌ DON'T: Use Fragile Selectors

```python
# Bad: CSS selectors (fragile)
page.locator(".btn-primary")
page.locator("#submit")

# Bad: XPath (fragile, hard to read)
page.locator("//button[@class='submit']")

# Bad: Text without specificity
page.get_by_text("Submit")  # May match multiple elements
```

### Handling Strict Mode Violations

When selectors match multiple elements, Playwright's strict mode raises an error.

```python
# Problem: "Retrosynthesis" matches 6 elements
page.get_by_text("Retrosynthesis Advisor")  # FAILS in strict mode

# Solution 1: Use .first (if first match is correct)
page.get_by_text("Retrosynthesis Advisor").first

# Solution 2: Scope to container
main = page.locator('[data-testid="stMainBlockContainer"]').first
main.get_by_text("Retrosynthesis Advisor")

# Solution 3: Use role with exact name (including emoji)
page.get_by_role("heading", name="🧪 Retrosynthesis Advisor")

# Solution 4: Chain locators
page.locator("main").get_by_role("button", name="Analyze")
```

## Streamlit-Specific Patterns

### Waiting for Streamlit Reruns

Streamlit reruns the entire script when inputs change. Use explicit waits:

```python
# After input, press Tab to trigger rerun
smiles_input.fill("CCO")
smiles_input.press("Tab")

# Wait for specific element to appear
page.get_by_role("button", name="Analyze").wait_for(state="visible")

# Wait for loading to complete
page.wait_for_load_state("domcontentloaded")
```

### Common Streamlit Locators

```python
# Main content area
main = page.locator('[data-testid="stMainBlockContainer"]').first

# Sidebar
sidebar = page.locator('[data-testid="stSidebar"]')

# Tabs
page.get_by_role("tab", name="Tab Name")

# Buttons
page.get_by_role("button", name="Button Text")

# Text inputs
page.get_by_label("Input Label")
# or
page.locator('input[aria-label="Input Label"]')
```

### Handling Dynamic Content

```python
# Wait for element to be visible
element = page.get_by_text("Results")
element.wait_for(state="visible", timeout=10000)

# Wait for element to disappear (loading spinner)
page.locator(".stSpinner").wait_for(state="hidden")

# Retry pattern for flaky elements
from playwright.sync_api import expect
expect(page.get_by_role("button", name="Submit")).to_be_enabled(timeout=5000)
```

## Test Structure

### Fixtures (conftest.py)

```python
import pytest
from playwright.sync_api import Page

@pytest.fixture
def base_url():
    """Streamlit server URL."""
    return "http://localhost:8501"

@pytest.fixture
def authenticated_page(page: Page, base_url: str):
    """Page with authentication bypassed."""
    # DISABLE_AUTH=true should be set in environment
    page.goto(base_url)
    page.wait_for_load_state("domcontentloaded")
    return page
```

### Test Class Pattern

```python
class TestFeatureE2E:
    """E2E tests for Feature dashboard."""

    def test_page_loads(self, page, base_url):
        """Dashboard page loads without errors."""
        page.goto(f"{base_url}/Feature")
        expect(page.get_by_role("heading", name="Feature")).to_be_visible()

    def test_form_submission(self, page, base_url):
        """Form submits and shows results."""
        page.goto(f"{base_url}/Feature")
        
        # Fill form
        page.get_by_label("Input").fill("test value")
        page.get_by_label("Input").press("Tab")
        
        # Submit
        page.get_by_role("button", name="Submit").click()
        
        # Verify results
        expect(page.get_by_text("Results")).to_be_visible()
```

## Debugging

### Headed Mode

```bash
# See browser during test
pytest test_file.py --headed

# Slow down for observation
pytest test_file.py --headed --slowmo=500
```

### Screenshots on Failure

```python
def test_something(page):
    try:
        # test code
    except Exception:
        page.screenshot(path="debug_screenshot.png")
        raise
```

### Trace Viewer

```bash
# Record trace
pytest test_file.py --tracing=on

# View trace
playwright show-trace trace.zip
```

## Common Issues

| Issue | Solution |
|-------|----------|
| Strict mode violation | Use `.first` or more specific selector |
| Element not found | Check if Streamlit rerun completed; add wait |
| Button disabled | Ensure all required inputs filled; press Tab |
| Timeout | Increase timeout or check if element exists |
| Auth redirect | Set `DISABLE_AUTH=true` in environment |

## Adding New E2E Tests

1. Create test file: `test_{feature}_e2e.py`
2. Use existing `conftest.py` fixtures
3. Follow patterns above for selectors
4. Run with `--headed` first to debug
5. Ensure tests pass in headless mode before commit

