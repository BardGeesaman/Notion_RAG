"""Quick actions utilities for command palette."""
from __future__ import annotations

from typing import List, Dict, Any

QUICK_ACTIONS = [
    {
        "name": "New Experiment",
        "description": "Create a new experiment",
        "page": "Experiments",
        "icon": "🧪",
    },
    {
        "name": "Register Compound",
        "description": "Register a new compound",
        "page": "Chemistry",
        "icon": "⚗️",
    },
    {
        "name": "Run Discovery",
        "description": "Run repository discovery job",
        "page": "Discovery Workflow",
        "icon": "🔍",
    },
    {
        "name": "Data Quality",
        "description": "Run data validation checks",
        "page": "Data Quality",
        "icon": "✅",
    },
    {
        "name": "Compare",
        "description": "Compare experiments or compounds",
        "page": "Compare",
        "icon": "🔍",
    },
    {
        "name": "Timeline",
        "description": "View activity timeline",
        "page": "Timeline",
        "icon": "📅",
    },
    {
        "name": "Import Data",
        "description": "Bulk import experiments, compounds, or samples",
        "page": "Import Data",
        "icon": "📥",
    },
    {
        "name": "Analysis Tools",
        "description": "Run pathway enrichment and analysis",
        "page": "Analysis Tools",
        "icon": "🔬",
    },
    {
        "name": "Search",
        "description": "Global search across entities",
        "page": "Search",
        "icon": "🔍",
    },
    {
        "name": "Overview",
        "description": "View platform overview and statistics",
        "page": "Overview",
        "icon": "📊",
    },
    {
        "name": "Protocols",
        "description": "Manage experimental protocols",
        "page": "Protocols",
        "icon": "📋",
    },
    {
        "name": "Sample Inventory",
        "description": "Manage sample inventory",
        "page": "Sample Inventory",
        "icon": "🧫",
    },
    {
        "name": "Q&A Tracker",
        "description": "Track questions and answers",
        "page": "Q&A Tracker",
        "icon": "❓",
    },
    {
        "name": "System Health",
        "description": "View system health and statistics",
        "page": "System Health",
        "icon": "🏥",
    },
]


def search_actions(query: str) -> List[Dict[str, Any]]:
    """
    Search quick actions by name or description (fuzzy match).

    Args:
        query: Search query string

    Returns:
        Filtered list of actions matching the query
    """
    if not query or not query.strip():
        return QUICK_ACTIONS

    query_lower = query.lower().strip()
    results = []

    for action in QUICK_ACTIONS:
        name_match = query_lower in action["name"].lower()
        desc_match = query_lower in action["description"].lower()

        if name_match or desc_match:
            results.append(action)

    return results
