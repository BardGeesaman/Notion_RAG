"""
Tests for Scientist's Cockpit dashboard widgets.

Tests widget rendering and API integration.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

# Import widget functions
from scripts.dashboard.components.cockpit_widgets import (
    render_stats_widget,
    render_activity_widget,
    render_alerts_widget,
    render_tasks_widget,
    render_shortcuts_widget,
)


class TestStatsWidget:
    """Tests for stats widget."""

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_stats_widget_success(self, mock_st, mock_api_get):
        """Test stats widget renders successfully."""
        mock_api_get.return_value = [{"id": "1"}]
        
        # Mock st.columns to return mock column objects
        mock_col1 = MagicMock()
        mock_col2 = MagicMock()
        mock_st.columns.return_value = [mock_col1, mock_col2]
        
        render_stats_widget()
        
        assert mock_st.markdown.called
        assert mock_st.columns.called

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_stats_widget_api_failure(self, mock_st, mock_api_get):
        """Test stats widget handles API failure."""
        mock_api_get.return_value = None
        
        mock_col1 = MagicMock()
        mock_col2 = MagicMock()
        mock_st.columns.return_value = [mock_col1, mock_col2]
        
        render_stats_widget()
        
        assert mock_st.markdown.called


class TestActivityWidget:
    """Tests for activity widget."""

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_activity_widget_with_events(self, mock_st, mock_api_get):
        """Test activity widget with events."""
        mock_api_get.return_value = [
            {"event_type": "experiment_created", "target_name": "Exp 1", "created_at": "2024-01-01"},
        ]
        
        render_activity_widget()
        
        assert mock_st.markdown.called

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_activity_widget_no_events(self, mock_st, mock_api_get):
        """Test activity widget with no events."""
        mock_api_get.return_value = []
        
        render_activity_widget()
        
        # Should show "No recent activity"
        assert mock_st.info.called or mock_st.markdown.called


class TestAlertsWidget:
    """Tests for alerts widget."""

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_alerts_widget_with_unread(self, mock_st, mock_api_get):
        """Test alerts widget with unread notifications."""
        mock_api_get.return_value = {"unread": 5}
        
        render_alerts_widget()
        
        assert mock_st.markdown.called
        # Should show warning for unread
        assert mock_st.warning.called or mock_st.success.called


class TestTasksWidget:
    """Tests for tasks widget."""

    @patch("scripts.dashboard.components.cockpit_widgets._api_get")
    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_tasks_widget_with_pending(self, mock_st, mock_api_get):
        """Test tasks widget with pending reviews."""
        mock_api_get.return_value = [
            {"notebook_path": "analysis.ipynb"},
        ]
        
        render_tasks_widget()
        
        assert mock_st.markdown.called


class TestShortcutsWidget:
    """Tests for shortcuts widget."""

    @patch("scripts.dashboard.components.cockpit_widgets.st")
    def test_render_shortcuts_widget(self, mock_st):
        """Test shortcuts widget renders buttons."""
        mock_col1 = MagicMock()
        mock_col2 = MagicMock()
        mock_st.columns.return_value = [mock_col1, mock_col2]
        
        render_shortcuts_widget()
        
        assert mock_st.markdown.called
        assert mock_st.columns.called

