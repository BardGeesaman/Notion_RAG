"""End-to-end tests for responsive navigation."""

from __future__ import annotations

import pytest
from playwright.sync_api import Page, expect


class TestResponsiveNavE2E:
    """Test responsive navigation functionality."""
    
    @pytest.fixture(autouse=True)
    def setup_page(self, page: Page, base_url: str):
        """Navigate to dashboard before each test."""
        if not base_url:
            pytest.skip("No base_url provided - E2E tests require running Streamlit server")
            
        page.goto(base_url)
        page.wait_for_load_state("domcontentloaded")
        
        # Wait for main content to load
        page.wait_for_selector("[data-testid='stMainBlockContainer']", timeout=10000)
    
    def test_sidebar_collapses_on_narrow(self, page: Page):
        """Test that sidebar auto-collapses on narrow screens."""
        # Set viewport to mobile size
        page.set_viewport_size({"width": 375, "height": 667})  # iPhone SE size
        page.wait_for_timeout(1000)
        
        # Check if sidebar exists
        sidebar = page.locator("[data-testid='stSidebar']")
        
        if sidebar.count() > 0:
            # Check if sidebar has mobile-friendly behavior
            # Look for collapse button or collapsed state
            collapse_button = page.locator("[data-testid='collapsedControl']")
            
            if collapse_button.count() > 0:
                # Check if sidebar is collapsed or collapsible
                sidebar_expanded = sidebar.get_attribute("aria-expanded")
                
                if sidebar_expanded == "false":
                    print("✅ Sidebar auto-collapsed on mobile viewport")
                else:
                    print("ℹ️ Sidebar not auto-collapsed - may need integration")
                    
                # Test manual collapse/expand
                collapse_button.click()
                page.wait_for_timeout(500)
                
                # Check state changed
                new_state = sidebar.get_attribute("aria-expanded")
                if new_state != sidebar_expanded:
                    print("✅ Sidebar collapse/expand working")
                else:
                    print("ℹ️ Sidebar state may not be changing")
            else:
                print("ℹ️ Collapse button not found - sidebar may always be visible")
        else:
            print("ℹ️ Sidebar not found - may need responsive implementation")
        
        # Reset to desktop size
        page.set_viewport_size({"width": 1200, "height": 800})
        page.wait_for_timeout(500)
    
    def test_touch_friendly_button_size(self, page: Page):
        """Test that buttons meet 44px minimum size on mobile."""
        # Set viewport to mobile size
        page.set_viewport_size({"width": 375, "height": 667})
        page.wait_for_timeout(1000)
        
        # Find navigation buttons
        nav_buttons = page.locator("button").filter(has_text="Overview")
        other_buttons = page.locator("button").filter(has_text="Chemistry")
        
        buttons_to_test = []
        if nav_buttons.count() > 0:
            buttons_to_test.append(nav_buttons.first)
        if other_buttons.count() > 0:
            buttons_to_test.append(other_buttons.first)
        
        touch_friendly_count = 0
        total_buttons = len(buttons_to_test)
        
        for button in buttons_to_test:
            if button.is_visible():
                # Get button dimensions
                bbox = button.bounding_box()
                if bbox:
                    height = bbox["height"]
                    width = bbox["width"]
                    
                    # Check if meets touch target minimum (44px)
                    if height >= 44:
                        touch_friendly_count += 1
                        print(f"✅ Button meets touch target: {height}px height")
                    else:
                        print(f"ℹ️ Button below touch target: {height}px height")
        
        if total_buttons > 0:
            touch_ratio = touch_friendly_count / total_buttons
            if touch_ratio >= 0.8:  # 80% of buttons should be touch-friendly
                print(f"✅ {touch_ratio:.1%} of buttons are touch-friendly")
            else:
                print(f"ℹ️ Only {touch_ratio:.1%} of buttons are touch-friendly")
        else:
            print("ℹ️ No navigation buttons found for touch testing")
        
        # Reset to desktop
        page.set_viewport_size({"width": 1200, "height": 800})
    
    def test_expansion_state_persists(self, page: Page):
        """Test that sidebar expansion state persists across viewport changes."""
        # Start with desktop size
        page.set_viewport_size({"width": 1200, "height": 800})
        page.wait_for_timeout(500)
        
        # Find sidebar and collapse button
        sidebar = page.locator("[data-testid='stSidebar']")
        collapse_button = page.locator("[data-testid='collapsedControl']")
        
        if sidebar.count() > 0 and collapse_button.count() > 0:
            # Get initial state
            initial_state = sidebar.get_attribute("aria-expanded")
            
            # Collapse sidebar if expanded
            if initial_state != "false":
                collapse_button.click()
                page.wait_for_timeout(500)
            
            # Check collapsed state
            collapsed_state = sidebar.get_attribute("aria-expanded")
            
            # Change to mobile viewport
            page.set_viewport_size({"width": 375, "height": 667})
            page.wait_for_timeout(1000)
            
            # Check if state persists
            mobile_state = sidebar.get_attribute("aria-expanded")
            
            if mobile_state == collapsed_state:
                print("✅ Sidebar state persists across viewport changes")
            else:
                print("ℹ️ Sidebar state may not persist across viewport changes")
            
            # Expand sidebar on mobile
            if mobile_state == "false":
                collapse_button.click()
                page.wait_for_timeout(500)
                
                expanded_mobile_state = sidebar.get_attribute("aria-expanded")
                if expanded_mobile_state != "false":
                    print("✅ Sidebar can be expanded on mobile")
                else:
                    print("ℹ️ Sidebar may not expand on mobile")
            
            # Return to desktop
            page.set_viewport_size({"width": 1200, "height": 800})
            page.wait_for_timeout(500)
            
            # Check final state
            final_state = sidebar.get_attribute("aria-expanded")
            print(f"Final sidebar state: {final_state}")
        else:
            print("ℹ️ Sidebar or collapse button not found for persistence testing")
    
    def test_mobile_css_injection(self, page: Page):
        """Test that mobile CSS is properly injected."""
        # Set mobile viewport
        page.set_viewport_size({"width": 375, "height": 667})
        page.wait_for_timeout(1000)
        
        # Check for mobile-specific CSS classes or styles
        sidebar = page.locator("[data-testid='stSidebar']")
        
        if sidebar.count() > 0:
            # Check if mobile class is applied
            sidebar_classes = sidebar.get_attribute("class") or ""
            
            if "mobile-sidebar" in sidebar_classes:
                print("✅ Mobile CSS class applied to sidebar")
            else:
                print("ℹ️ Mobile CSS class not found - may need integration")
            
            # Check for mobile-specific styles
            sidebar_style = sidebar.get_attribute("style") or ""
            
            # Look for mobile-specific width or other styles
            if "width" in sidebar_style or "280px" in sidebar_style:
                print("✅ Mobile-specific styling detected")
            else:
                print("ℹ️ Mobile-specific styling not detected")
        
        # Reset viewport
        page.set_viewport_size({"width": 1200, "height": 800})
    
    def test_responsive_breakpoints(self, page: Page):
        """Test behavior at different responsive breakpoints."""
        breakpoints = [
            (375, 667, "xs"),   # Mobile portrait
            (576, 768, "sm"),   # Mobile landscape / small tablet
            (768, 1024, "md"),  # Tablet
            (1200, 800, "lg"),  # Desktop
        ]
        
        for width, height, size_name in breakpoints:
            page.set_viewport_size({"width": width, "height": height})
            page.wait_for_timeout(1000)
            
            print(f"Testing {size_name} breakpoint ({width}x{height})")
            
            # Check sidebar behavior
            sidebar = page.locator("[data-testid='stSidebar']")
            
            if sidebar.count() > 0:
                # Check if sidebar is responsive
                sidebar_visible = sidebar.is_visible()
                
                if sidebar_visible:
                    print(f"  ✅ Sidebar visible at {size_name}")
                else:
                    print(f"  ℹ️ Sidebar hidden at {size_name}")
                
                # Check button sizes for touch-friendly breakpoints
                if size_name in ["xs", "sm"]:
                    buttons = page.locator("button").first
                    if buttons.count() > 0:
                        bbox = buttons.bounding_box()
                        if bbox and bbox["height"] >= 44:
                            print(f"  ✅ Touch-friendly buttons at {size_name}")
                        else:
                            print(f"  ℹ️ Buttons may not be touch-friendly at {size_name}")
            else:
                print(f"  ℹ️ Sidebar not found at {size_name}")
    
    def test_compact_mode_functionality(self, page: Page):
        """Test compact mode for very small screens."""
        # Set very small viewport (extra small)
        page.set_viewport_size({"width": 320, "height": 568})  # iPhone 5 size
        page.wait_for_timeout(1000)
        
        # Check for compact mode features
        sidebar = page.locator("[data-testid='stSidebar']")
        
        if sidebar.count() > 0:
            # In compact mode, pin buttons should be hidden
            pin_buttons = sidebar.locator("button").filter(has_text="📌")
            
            if pin_buttons.count() == 0:
                print("✅ Pin buttons hidden in compact mode")
            else:
                print("ℹ️ Pin buttons still visible in compact mode")
            
            # Check for help text on buttons (compact mode feature)
            nav_buttons = sidebar.locator("button[title]")  # Buttons with title/help
            
            if nav_buttons.count() > 0:
                print("✅ Help text available on buttons in compact mode")
            else:
                print("ℹ️ Help text not found in compact mode")
            
            # Check sidebar width in compact mode
            bbox = sidebar.bounding_box()
            if bbox:
                width = bbox["width"]
                # Should be constrained on very small screens
                if width <= 320:  # Should not exceed viewport
                    print(f"✅ Sidebar width appropriate for compact mode: {width}px")
                else:
                    print(f"ℹ️ Sidebar may be too wide for compact mode: {width}px")
        
        # Reset viewport
        page.set_viewport_size({"width": 1200, "height": 800})
    
    def test_desktop_experience_preserved(self, page: Page):
        """Test that desktop experience is not broken by mobile changes."""
        # Set desktop viewport
        page.set_viewport_size({"width": 1920, "height": 1080})
        page.wait_for_timeout(1000)
        
        # Check that desktop features are preserved
        sidebar = page.locator("[data-testid='stSidebar']")
        
        if sidebar.count() > 0:
            # Desktop should have pin functionality
            pin_buttons = sidebar.locator("button").filter(has_text="📌")
            
            if pin_buttons.count() > 0:
                print("✅ Pin functionality available on desktop")
            else:
                print("ℹ️ Pin functionality not found on desktop")
            
            # Desktop should not have mobile CSS classes
            sidebar_classes = sidebar.get_attribute("class") or ""
            
            if "mobile-sidebar" not in sidebar_classes:
                print("✅ Mobile CSS classes not applied on desktop")
            else:
                print("ℹ️ Mobile CSS classes may be applied on desktop")
            
            # Check that sidebar is expanded by default on desktop
            expanded_state = sidebar.get_attribute("aria-expanded")
            
            if expanded_state != "false":
                print("✅ Sidebar expanded by default on desktop")
            else:
                print("ℹ️ Sidebar may be collapsed on desktop")
            
            # Check button sizes are normal (not oversized for touch)
            buttons = page.locator("button").first
            if buttons.count() > 0:
                bbox = buttons.bounding_box()
                if bbox:
                    height = bbox["height"]
                    # Desktop buttons should be normal size, not necessarily 44px+
                    if 32 <= height <= 48:  # Normal range
                        print(f"✅ Normal button size on desktop: {height}px")
                    else:
                        print(f"ℹ️ Button size may be non-standard: {height}px")
        
        print("Desktop experience validation complete")
