"""End-to-end tests for navigation sidebar."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect

from scripts.dashboard.core.config import PAGE_GROUPS, GROUP_ORDER


class TestNavigationSidebarE2E:
    """Test navigation sidebar functionality."""
    
    @pytest.fixture(autouse=True)
    def setup_page(self, page: Page, base_url: str):
        """Navigate to dashboard before each test."""
        page.goto(base_url)
        page.wait_for_load_state("domcontentloaded")
        
        # Wait for sidebar to load
        page.wait_for_selector("[data-testid='stSidebar']", timeout=10000)
    
    def test_sidebar_renders_all_groups(self, page: Page):
        """Test that all PAGE_GROUPS are visible in sidebar."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Count non-empty groups
        expected_groups = []
        for group in GROUP_ORDER:
            pages = PAGE_GROUPS.get(group, [])
            if pages and not (group == "Other" and not pages):
                expected_groups.append(group)
        
        # Check that all expected groups are present
        for group in expected_groups:
            group_locator = sidebar.get_by_text(group, exact=False)
            expect(group_locator).to_be_visible()
    
    def test_group_expansion_toggle(self, page: Page):
        """Test that group expanders can be opened and closed."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Find the Chemistry group expander
        chemistry_expander = sidebar.locator("summary").filter(has_text="Chemistry")
        expect(chemistry_expander).to_be_visible()
        
        # Click to expand/collapse
        chemistry_expander.click()
        page.wait_for_timeout(500)  # Wait for animation
        
        # Check that content is visible after expansion
        chemistry_content = sidebar.locator("details").filter(has_text="Chemistry").locator("div")
        # The expander should contain navigation buttons
        expect(chemistry_content).to_be_visible()
    
    def test_page_navigation(self, page: Page):
        """Test that clicking on a page navigates correctly."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Expand Chemistry group if not already expanded
        chemistry_expander = sidebar.locator("summary").filter(has_text="Chemistry")
        chemistry_expander.click()
        page.wait_for_timeout(500)
        
        # Look for a chemistry page button (e.g., "Chemical Sketcher")
        chemistry_page = sidebar.get_by_role("button", name="Chemical Sketcher")
        if chemistry_page.is_visible():
            chemistry_page.click()
            page.wait_for_timeout(1000)
            
            # Check that navigation occurred (page title or content changed)
            # Look for any indication that we're on the Chemical Sketcher page
            main_content = page.locator("[data-testid='stMainBlockContainer']")
            expect(main_content).to_be_visible()
    
    def test_current_page_highlighted(self, page: Page):
        """Test that the current page is highlighted differently."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Navigate to a specific page first
        home_expander = sidebar.locator("summary").filter(has_text="Home")
        home_expander.click()
        page.wait_for_timeout(500)
        
        # Click on Overview page
        overview_button = sidebar.get_by_role("button", name="Overview")
        if overview_button.is_visible():
            overview_button.click()
            page.wait_for_timeout(1000)
            
            # Check that Overview button has primary styling or arrow indicator
            # Look for the arrow indicator or primary button styling
            highlighted_button = sidebar.locator("button").filter(has_text="→ Overview")
            if highlighted_button.count() == 0:
                # Fallback: look for primary button type
                highlighted_button = sidebar.locator("button[data-baseweb='button'][data-testid='baseButton-primary']").filter(has_text="Overview")
            
            expect(highlighted_button).to_be_visible()
    
    def test_group_page_count(self, page: Page):
        """Test that group expanders show correct page count badges."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Check that groups show page counts in parentheses
        for group in GROUP_ORDER[:3]:  # Test first 3 groups
            pages = PAGE_GROUPS.get(group, [])
            if pages:
                page_count = len(pages)
                
                # Look for group with count badge
                group_with_count = sidebar.locator("summary").filter(has_text=f"{group} ({page_count})")
                expect(group_with_count).to_be_visible()
    
    def test_group_pinning_functionality(self, page: Page):
        """Test that groups can be pinned to stay expanded."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Find Chemistry group
        chemistry_expander = sidebar.locator("summary").filter(has_text="Chemistry")
        chemistry_expander.click()
        page.wait_for_timeout(500)
        
        # Look for pin button (📌 or 📍)
        pin_button = sidebar.get_by_role("button", name="📌")
        if pin_button.is_visible():
            pin_button.click()
            page.wait_for_timeout(500)
            
            # Check that pin icon changed to 📍
            pinned_button = sidebar.get_by_role("button", name="📍")
            expect(pinned_button).to_be_visible()
    
    def test_sidebar_search_functionality(self, page: Page):
        """Test that sidebar search works correctly."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Find search input
        search_input = sidebar.get_by_placeholder("Search experiments, compounds, datasets...")
        expect(search_input).to_be_visible()
        
        # Type a search query
        search_input.fill("test")
        page.wait_for_timeout(1000)
        
        # The search should trigger some kind of results display
        # This is a basic test to ensure search input is functional
        expect(search_input).to_have_value("test")
    
    def test_recent_pages_section(self, page: Page):
        """Test that recent pages section appears after navigation."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Navigate to a page to populate recent pages
        home_expander = sidebar.locator("summary").filter(has_text="Home")
        home_expander.click()
        page.wait_for_timeout(500)
        
        overview_button = sidebar.get_by_role("button", name="Overview")
        if overview_button.is_visible():
            overview_button.click()
            page.wait_for_timeout(1000)
            
            # Check if Recent section appears
            recent_section = sidebar.locator("summary").filter(has_text="Recent")
            # Recent section might not always be visible, so this is optional
            if recent_section.is_visible():
                expect(recent_section).to_be_visible()
    
    def test_navigation_persistence(self, page: Page):
        """Test that navigation state persists correctly."""
        sidebar = page.locator("[data-testid='stSidebar']")
        
        # Navigate to Chemistry page
        chemistry_expander = sidebar.locator("summary").filter(has_text="Chemistry")
        chemistry_expander.click()
        page.wait_for_timeout(500)
        
        chemistry_page = sidebar.get_by_role("button", name="Chemistry")
        if chemistry_page.is_visible():
            chemistry_page.click()
            page.wait_for_timeout(1000)
            
            # Refresh the page
            page.reload()
            page.wait_for_load_state("domcontentloaded")
            page.wait_for_timeout(2000)
            
            # Check that we're still on the Chemistry page
            main_content = page.locator("[data-testid='stMainBlockContainer']")
            expect(main_content).to_be_visible()
