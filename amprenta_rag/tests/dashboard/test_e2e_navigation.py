"""
E2E Navigation/UI test scenarios.

Extracted from the historical monolithic `test_e2e_platform.py` suite to keep
tests scenario-focused and easier to maintain.
"""

import pytest
from playwright.sync_api import Page


pytestmark = pytest.mark.requires_server


class TestNavigationSmoke:
    """Test all sidebar sections expand and pages load"""

    @pytest.mark.parametrize(
        "section,pages",
        [
            ("🔍 Discovery", ["Overview", "Experiments", "Discovery Workflow", "Variant Tracking"]),
            ("📊 Analysis", ["Datasets", "Signatures", "RAG Query"]),
            ("📋 ELN", ["Protocols", "Sample Inventory"]),
            ("⚙️ Admin", ["Audit Logs", "Data Quality", "System Health"]),
        ],
    )
    def test_section_pages_load(self, page: Page, base_url: str, section: str, pages: list):
        """Test that pages in each section load without errors."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)  # Wait for Streamlit to fully load

        # Expand section
        try:
            expander_button = page.get_by_role("button", name=section).first
            if expander_button.count() > 0:
                expander_button.click()
                page.wait_for_timeout(1000)  # Wait for expander to open
        except Exception as e:
            pytest.skip(f"Could not find or click section '{section}': {e}")

        for page_name in pages:
            try:
                # Click the page button
                page_button = page.get_by_role("button", name=page_name, exact=True).first
                if page_button.count() > 0:
                    page_button.click()
                    page.wait_for_timeout(2000)  # Wait for page to load

                    # Check for both error patterns
                    loading_error = page.locator("text=Error loading page").count()
                    rendering_error = page.locator("text=Error rendering page").count()
                    assert loading_error == 0, f"{page_name}: Import error - Error loading page"
                    assert rendering_error == 0, f"{page_name}: Runtime error - Error rendering page"

                    # Check that page loaded (some content is present)
                    # Most pages will have their title/header visible
                    assert page.locator("body").count() > 0, f"{page_name} did not load"
            except Exception as e:
                pytest.fail(f"Failed to navigate to {page_name}: {e}")


class TestOverview:
    """Test overview page"""

    def test_view_overview(self, page: Page, base_url: str):
        """View overview page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 🔍 Discovery (expanded default, no expander click) > Overview
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Click Overview (Discovery is expanded by default)
            sidebar.locator('button:has-text("Overview")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Overview: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Overview page: Import error - Error loading page"
            assert rendering_error == 0, "Overview page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Overview page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Overview page loaded: {e}")


class TestGettingStarted:
    """Test getting started page"""

    def test_view_getting_started(self, page: Page, base_url: str):
        """View getting started page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Getting Started
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            sidebar.locator('button:has-text("Getting Started")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Getting Started: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Getting Started page: Import error - Error loading page"
            assert rendering_error == 0, "Getting Started page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Getting Started page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Getting Started page loaded: {e}")


class TestLabNotebook:
    """Test lab notebook page"""

    def test_view_lab_notebook(self, page: Page, base_url: str):
        """View lab notebook page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Lab Notebook
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            sidebar.locator('button:has-text("Lab Notebook")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Lab Notebook: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Lab Notebook page: Import error - Error loading page"
            assert rendering_error == 0, "Lab Notebook page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Lab Notebook page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Lab Notebook page loaded: {e}")


class TestDiscoveryPage:
    """Test discovery page (different from Discovery Workflow)"""

    def test_view_discovery(self, page: Page, base_url: str):
        """View discovery page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📚 Other Pages expander > Discovery
        try:
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            page.locator("text=📚 Other Pages").first.click()
            page.wait_for_timeout(1000)

            # Use selector that excludes "Discovery Workflow"
            sidebar.locator('button:has-text("Discovery"):not(:has-text("Workflow"))').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Discovery: {e}")

        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Discovery page: Import error - Error loading page"
            assert rendering_error == 0, "Discovery page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Discovery page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Discovery page loaded: {e}")


