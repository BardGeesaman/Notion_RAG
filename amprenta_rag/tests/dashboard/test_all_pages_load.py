from __future__ import annotations

"""
Legacy AppTest-based page rendering tests have been removed due to
collection-time monkeypatch issues. Page-level behavior is now covered by:
- Import tests in test_pages_import.py
- Playwright-based E2E/UI tests
"""

import pytest


@pytest.mark.ui
def test_app_test_replaced_by_import_and_playwright():
    pytest.skip(
        "AppTest-based rendering has been replaced by import-only and Playwright UI tests."
    )


