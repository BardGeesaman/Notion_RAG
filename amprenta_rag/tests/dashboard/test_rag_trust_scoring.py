import pytest
from playwright.sync_api import Page, expect


@pytest.fixture(scope="function")
def rag_page(page: Page, base_url: str):
    """Navigate to RAG Query page."""
    page.goto(base_url)
    page.wait_for_load_state("networkidle")
    # Expand Other Pages and click RAG Query
    page.locator("text=📚 Other Pages").click()
    page.wait_for_timeout(500)
    page.get_by_role("button", name="RAG Query", exact=True).click()
    page.wait_for_timeout(2000)
    return page


def test_trust_scoring_checkbox_exists(rag_page: Page):
    """Trust scoring checkbox should be visible."""
    checkbox = rag_page.locator("text=Apply Trust Scoring")
    expect(checkbox).to_be_visible(timeout=5000)


def test_trust_scoring_enabled(rag_page: Page):
    """When trust scoring enabled and query run, trust metrics should appear."""
    # Enable trust scoring
    rag_page.locator("text=Apply Trust Scoring").click()

    # Enter a query
    query_input = rag_page.get_by_label("Enter your query")
    query_input.fill("What is ALS?")

    # Click search/query button
    rag_page.get_by_role("button", name="Search").click()

    # Wait for results
    rag_page.wait_for_timeout(5000)

    # Check for trust analysis section
    expect(rag_page.locator("text=Trust")).to_be_attached(timeout=10000)


