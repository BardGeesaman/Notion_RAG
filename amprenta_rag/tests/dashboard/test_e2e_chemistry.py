"""
E2E Chemistry test scenarios.

Extracted from the historical monolithic `test_e2e_platform.py` suite.
"""

import time

import pytest
from playwright.sync_api import Page, expect


pytestmark = pytest.mark.requires_server


class TestChemistryFlow:
    """Test compound registration and viewing"""

    def test_register_and_view_compound(self, page: Page, base_url: str):
        """Test registering a compound and verifying it appears in the list."""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate to Chemistry
        try:
            # Scroll sidebar to reveal expanders
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 400)
            page.wait_for_timeout(1000)

            # Now click Analysis expander (Streamlit expander, not a button)
            page.locator("text=📊 Analysis").click()
            page.wait_for_timeout(1000)
            # Then click Chemistry
            page.get_by_role("button", name="Chemistry", exact=True).click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Chemistry: {e}")

        # Go to Register tab (use text locator for Streamlit tabs)
        try:
            page.locator('[data-testid="stTabs"] button:has-text("Register Compound")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not find Register Compound tab: {e}")

        # Fill form with unique compound name
        compound_name = f"TestCompound_E2E_{int(time.time())}"

        try:
            # Use exact aria-label attribute selectors to avoid partial matches
            name_input = page.locator('input[aria-label="Compound Name"]')
            name_input.fill(compound_name)

            smiles_input = page.locator('input[aria-label="SMILES"]')
            smiles_input.fill("CCO")  # Ethanol

            page.wait_for_timeout(500)

            # Click Register button
            register_button = page.get_by_role("button", name="Register", exact=True).first
            if register_button.count() > 0:
                register_button.click()
            else:
                # Try without exact match
                page.get_by_role("button", name="Register").first.click()

            page.wait_for_timeout(3000)  # Wait for registration to complete
        except Exception as e:
            pytest.fail(f"Failed to fill form or register compound: {e}")

        # Verify registration (check for AMP- ID pattern)
        try:
            amp_id_locator = page.locator("text=AMP-").first
            expect(amp_id_locator).to_be_attached(timeout=10000)
        except Exception as e:
            pytest.fail(f"Compound registration verification failed: {e}")


class TestCandidateSelection:
    """Test candidate selection page"""

    def test_view_candidate_selection(self, page: Page, base_url: str):
        """View candidate selection page"""
        page.goto(base_url)
        page.wait_for_load_state("networkidle")
        page.wait_for_timeout(3000)

        # Navigate: 📊 Analysis expander > Candidate Selection
        try:
            # Sidebar setup
            sidebar = page.locator('[data-testid="stSidebar"]').first
            sidebar.hover()
            page.mouse.wheel(0, 600)
            page.wait_for_timeout(2000)

            # Expand Analysis section
            page.locator("text=📊 Analysis").first.click()
            page.wait_for_timeout(1000)

            # Click Candidate Selection
            sidebar.locator('button:has-text("Candidate Selection")').click()
            page.wait_for_timeout(2000)
        except Exception as e:
            pytest.fail(f"Could not navigate to Candidate Selection: {e}")

        # Verify page loads without error
        try:
            loading_error = page.locator("text=Error loading page").count()
            rendering_error = page.locator("text=Error rendering page").count()
            assert loading_error == 0, "Candidate Selection page: Import error - Error loading page"
            assert rendering_error == 0, "Candidate Selection page: Runtime error - Error rendering page"
            assert page.locator("body").count() > 0, "Candidate Selection page did not load"
        except Exception as e:
            pytest.fail(f"Could not verify Candidate Selection page loaded: {e}")


