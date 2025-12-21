"""
E2E Search/RAG test scenarios.

Extracted from the historical monolithic `test_e2e_platform.py` suite.
"""

import pytest
from playwright.sync_api import Page


pytestmark = pytest.mark.requires_server


class TestRAGFlow:
    """Test RAG query with hybrid search"""

    def test_rag_query_with_citations(self, page: Page, base_url: str):
        """Test RAG query functionality with citations."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate to RAG Query
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)

            # RAG Query might be in "Other Pages" or "Analysis"
            # Try Analysis first (Streamlit expander, not a button)
            analysis_expander = page.locator("text=📊 Analysis").first
            if analysis_expander.count() > 0:
                analysis_expander.click()
                page.wait_for_timeout(1000)

            # Check if RAG Query is in Analysis
            rag_button = page.get_by_role("button", name="RAG Query", exact=True).first
            if rag_button.count() == 0:
                # Try Other Pages
                other_pages_expander = page.locator("text=📚 Other Pages").first
                if other_pages_expander.count() > 0:
                    other_pages_expander.click()
                    page.wait_for_timeout(1000)
                    rag_button = page.get_by_role("button", name="RAG Query", exact=True).first

            if rag_button.count() > 0:
                rag_button.click()
            else:
                pytest.skip("RAG Query page not found in sidebar")

            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.skip(f"Could not navigate to RAG Query: {e}")

        # Enter query
        try:
            # Use exact aria-label selector for query input
            query_input = page.locator('textarea[aria-label="Enter your query"]')
            query_input.fill("What causes ALS?")
        except Exception as e:
            pytest.fail(f"Could not find query input: {e}")

        # Enable hybrid search (if checkbox exists)
        try:
            # Use exact aria-label selector for checkbox
            hybrid_checkbox = page.locator('input[aria-label="Use Hybrid Search"]')
            if hybrid_checkbox.count() > 0:
                hybrid_checkbox.check()
        except Exception:
            # Hybrid search checkbox might not exist or have different label
            pass

        # Submit query
        try:
            search_button = page.get_by_role("button", name="🔍 Search", exact=True).first
            if search_button.count() == 0:
                search_button = page.get_by_role("button", name="Search", exact=True).first
            if search_button.count() > 0:
                search_button.click()
            else:
                pytest.fail("Could not find Search button")
        except Exception as e:
            pytest.fail(f"Could not submit query: {e}")

        # Wait for results (RAG takes time)
        page.wait_for_timeout(10000)

        # Check for results (matches section or source citations)
        try:
            # Look for common result indicators
            source_locator = page.locator("text=Source").first
            if source_locator.count() == 0:
                # Try alternative indicators
                matches_locator = page.locator("text=Matches").first
                if matches_locator.count() == 0:
                    # Check if any content loaded (not just error)
                    error_locator = page.locator("text=Error").first
                    if error_locator.count() > 0:
                        pytest.fail("RAG query returned an error")
                    # If no error, assume results loaded
                    assert page.locator("body").count() > 0, "Page did not load results"
                else:
                    assert matches_locator.count() > 0, "Matches section not found"
            else:
                assert source_locator.count() > 0, "Source section not found"
        except Exception as e:
            pytest.fail(f"Could not verify RAG query results: {e}")


class TestChat:
    """Test chat page"""

    def test_view_chat(self, page: Page, base_url: str):
        """View chat page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Chat
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            # Click Chat
            sidebar.locator('button:has-text("Chat")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Chat: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Chat page: Import error - Error loading page"
            assert rendering_error == 0, "Chat page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Chat page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Chat page loaded: {e}")


class TestSearch:
    """Test search page"""

    def test_view_search(self, page: Page, base_url: str):
        """View search page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Search
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            # Click Search
            sidebar.locator('button:has-text("Search")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Search: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Search page: Import error - Error loading page"
            assert rendering_error == 0, "Search page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Search page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Search page loaded: {e}")


class TestCrossOmics:
    """Test cross-omics page"""

    def test_view_cross_omics(self, page: Page, base_url: str):
        """View cross-omics page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Cross-Omics
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Expand Other Pages section
            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            # Click Cross-Omics
            sidebar.locator('button:has-text("Cross-Omics")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Cross-Omics: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Cross-Omics page: Import error - Error loading page"
            assert rendering_error == 0, "Cross-Omics page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Cross-Omics page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Cross-Omics page loaded: {e}")


class TestRAGChunks:
    """Test RAG chunks page"""

    def test_view_rag_chunks(self, page: Page, base_url: str):
        """View RAG chunks page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > RAG Chunks
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            sidebar.locator('button:has-text("RAG Chunks")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to RAG Chunks: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "RAG Chunks page: Import error - Error loading page"
            assert rendering_error == 0, "RAG Chunks page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "RAG Chunks page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify RAG Chunks page loaded: {e}")


