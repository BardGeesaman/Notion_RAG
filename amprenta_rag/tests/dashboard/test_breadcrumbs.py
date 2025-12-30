"""Tests for breadcrumb navigation component."""

from __future__ import annotations

import pytest

from scripts.dashboard.components.breadcrumbs import (
    get_breadcrumb_path,
    get_page_group,
    get_group_icon,
)
from scripts.dashboard.core.config import PAGE_GROUPS, GROUP_ICONS


class TestBreadcrumbs:
    """Test breadcrumb navigation functionality."""
    
    def test_breadcrumb_path_chemistry_page(self):
        """Test breadcrumb path for a Chemistry group page."""
        path = get_breadcrumb_path("Generative Chemistry")
        
        expected = [
            ("🏠", "Home", "Overview"),
            ("⚗️", "Chemistry", None),
            ("", "Generative Chemistry", "Generative Chemistry"),
        ]
        
        assert path == expected
    
    def test_breadcrumb_path_admin_page(self):
        """Test breadcrumb path for an Admin group page."""
        path = get_breadcrumb_path("Company Settings")
        
        expected = [
            ("🏠", "Home", "Overview"),
            ("⚙️", "Admin", None),
            ("", "Company Settings", "Company Settings"),
        ]
        
        assert path == expected
    
    def test_breadcrumb_path_home_page(self):
        """Test breadcrumb path for Home group pages."""
        # Test Overview page (special case)
        overview_path = get_breadcrumb_path("Overview")
        expected_overview = [
            ("🏠", "Home", None),
            ("", "Overview", "Overview"),
        ]
        assert overview_path == expected_overview
        
        # Test other Home page
        cockpit_path = get_breadcrumb_path("Cockpit")
        expected_cockpit = [
            ("🏠", "Home", "Overview"),
            ("", "Cockpit", "Cockpit"),
        ]
        assert cockpit_path == expected_cockpit
    
    def test_breadcrumb_path_unknown_page(self):
        """Test breadcrumb path fallback for pages not in groups."""
        path = get_breadcrumb_path("Unknown Page")
        
        expected = [
            ("🏠", "Home", "Overview"),
            ("", "Unknown Page", "Unknown Page"),
        ]
        
        assert path == expected
    
    def test_breadcrumb_path_discovery_page(self):
        """Test breadcrumb path for a Discovery group page."""
        path = get_breadcrumb_path("Paper Search")
        
        expected = [
            ("🏠", "Home", "Overview"),
            ("🔍", "Discovery", None),
            ("", "Paper Search", "Paper Search"),
        ]
        
        assert path == expected
    
    def test_breadcrumb_path_structure(self):
        """Test that breadcrumb paths have correct structure."""
        test_pages = ["Chemistry", "Audit Logs", "Overview", "Notebooks"]
        
        for page in test_pages:
            path = get_breadcrumb_path(page)
            
            # All paths should be lists of tuples
            assert isinstance(path, list)
            
            for item in path:
                assert isinstance(item, tuple)
                assert len(item) == 3
                
                icon, label, page_name = item
                assert isinstance(icon, str)
                assert isinstance(label, str)
                assert page_name is None or isinstance(page_name, str)
    
    def test_get_page_group_function(self):
        """Test the get_page_group helper function."""
        # Test known pages
        assert get_page_group("Chemical Sketcher") == "Chemistry"
        assert get_page_group("Company Settings") == "Admin"
        assert get_page_group("Overview") == "Home"
        assert get_page_group("Paper Search") == "Discovery"
        
        # Test unknown page
        assert get_page_group("Unknown Page") is None
    
    def test_get_group_icon_function(self):
        """Test the get_group_icon helper function."""
        # Test known groups
        assert get_group_icon("Chemistry") == "⚗️"
        assert get_group_icon("Admin") == "⚙️"
        assert get_group_icon("Home") == "🏠"
        assert get_group_icon("Discovery") == "🔍"
        
        # Test unknown group (should return default)
        assert get_group_icon("Unknown Group") == "📁"
    
    def test_breadcrumb_path_covers_all_groups(self):
        """Test that breadcrumb paths work for all defined groups."""
        for group, pages in PAGE_GROUPS.items():
            if pages:  # Skip empty groups
                # Test first page in each group
                test_page = pages[0]
                path = get_breadcrumb_path(test_page)
                
                # Should have at least 2 items (Home + current page)
                assert len(path) >= 2
                
                # First item should always be Home
                assert path[0][1] == "Home"
                
                # If not a Home page, should have group in path
                if group != "Home":
                    assert len(path) == 3  # Home -> Group -> Page
                    assert path[1][1] == group
                    assert path[2][1] == test_page
    
    def test_breadcrumb_icons_consistency(self):
        """Test that breadcrumb icons match GROUP_ICONS configuration."""
        for group in PAGE_GROUPS.keys():
            if group in GROUP_ICONS:
                icon_from_config = GROUP_ICONS[group]
                icon_from_function = get_group_icon(group)
                assert icon_from_config == icon_from_function
    
    def test_breadcrumb_navigation_logic(self):
        """Test the navigation logic in breadcrumb paths."""
        path = get_breadcrumb_path("Generative Chemistry")
        
        # Home should be clickable (has page_name)
        home_item = path[0]
        assert home_item[2] == "Overview"  # Should navigate to Overview
        
        # Group should not be clickable (page_name is None)
        group_item = path[1]
        assert group_item[2] is None
        
        # Current page should have itself as page_name but won't be clickable in UI
        current_item = path[2]
        assert current_item[2] == "Generative Chemistry"
    
    def test_empty_page_name_handling(self):
        """Test handling of empty or None page names."""
        # Test with empty string
        path = get_breadcrumb_path("")
        expected = [
            ("🏠", "Home", "Overview"),
            ("", "", ""),
        ]
        assert path == expected
        
        # Function should handle gracefully without crashing
        assert isinstance(path, list)
        assert len(path) >= 1
