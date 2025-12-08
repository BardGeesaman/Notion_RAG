from __future__ import annotations

"""
Simple import tests for all Streamlit dashboard pages.

These tests ensure that each page module can be imported without raising
ImportError or other exceptions. Detailed UI behavior is covered by the
Playwright-based tests.
"""

import pytest


PAGES = [
    "overview",
    "getting_started",
    "evaluation_wizard",
    "chat",
    "lab_notebook",
    "search",
    "ingestion",
    "repositories",
    "analysis",
    "discovery",
    "coverage",
    "feature_recurrence",
    "evidence_report",
    "management",
    "health",
    "relationships",
    "datasets",
    "programs",
    "experiments",
    "features",
    "signatures",
    "literature",
    "emails",
    "rag_chunks",
    "chemistry",
    "rag_query",
    "cross_omics",
]


@pytest.mark.ui
@pytest.mark.parametrize("page_name", PAGES)
def test_page_imports(page_name: str):
    module_name = f"scripts.dashboard.pages.{page_name}"
    __import__(module_name)


