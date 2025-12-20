"""Help content for dashboard pages."""
from __future__ import annotations

from typing import Dict, List, Optional, Any

HELP_TOPICS: Dict[str, Dict[str, Any]] = {
    "Overview": {
        "title": "Overview Dashboard",
        "description": "View system statistics, recent activity, and quick access to key features.",
        "tips": [
            "Use the metrics cards to quickly see dataset and experiment counts",
            "Recent experiments table shows latest additions to the system",
            "Click on any experiment name to view details",
        ],
    },
    "Experiments": {
        "title": "Experiments Browser",
        "description": "Browse and search all experiments in the system. Filter by disease, design type, or metadata.",
        "tips": [
            "Use the search bar to find experiments by name or description",
            "Filter by design type (case_control, time_course, etc.) to narrow results",
            "Click 'View Details' to see full experiment metadata and linked datasets",
        ],
    },
    "Chemistry": {
        "title": "Chemistry Informatics",
        "description": "Manage compounds, register new molecules, perform structure searches, and analyze SAR data.",
        "tips": [
            "Register compounds with SMILES strings to get automatic corporate IDs",
            "Use substructure search to find compounds matching a SMARTS pattern",
            "SAR Analysis tab shows property distributions and Lipinski Rule of 5 compliance",
            "Lead Optimization tab tracks ADME/PK/Tox data for compounds",
        ],
    },
    "Analysis Tools": {
        "title": "Analysis Tools",
        "description": "Perform pathway enrichment, signature scoring, and other analytical workflows.",
        "tips": [
            "Upload gene lists for pathway enrichment analysis",
            "Score signatures against datasets to find matches",
            "Export results as CSV for further analysis",
        ],
    },
    "Protocols": {
        "title": "Protocol Management",
        "description": "Create and manage experimental protocols. Link protocols to experiments for reproducibility.",
        "tips": [
            "Create reusable protocol templates with step-by-step instructions",
            "Link protocols to experiments to track which methods were used",
            "Version protocols by creating new versions from existing ones",
        ],
    },
    "Sample Inventory": {
        "title": "Sample Inventory",
        "description": "Track biological samples, storage locations, and sample transfers.",
        "tips": [
            "Create hierarchical storage locations (freezer → shelf → box)",
            "Register samples with barcodes for easy tracking",
            "Use sample transfers to record location changes",
            "View sample lineage to track parent-child relationships",
        ],
    },
    "Discovery Workflow": {
        "title": "Discovery Workflow",
        "description": "Automated repository scanning and study import. Schedule regular harvests from GEO and MetaboLights.",
        "tips": [
            "Run discovery jobs to scan repositories for new studies",
            "Review pending studies before importing as experiments",
            "Create harvest schedules for automated weekly/monthly scans",
            "Check job history to see past discovery runs",
        ],
    },
    "Q&A Tracker": {
        "title": "Q&A Tracker",
        "description": "Save questions and answers from RAG queries. Track question versions and export for documentation.",
        "tips": [
            "Ask questions and save both question and answer for future reference",
            "Re-run saved questions to get updated answers",
            "Export Q&A pairs as CSV for documentation",
            "Search your saved questions by keywords",
        ],
    },
    "Teams & Projects": {
        "title": "Teams & Projects",
        "description": "Organize work into teams and projects. Control access and collaboration.",
        "tips": [
            "Create teams to group collaborators",
            "Add team members with roles (owner, admin, member, viewer)",
            "Create projects within teams to organize experiments",
            "Make projects public to share across teams",
        ],
    },
    "Feedback": {
        "title": "Feedback & Feature Requests",
        "description": "Submit bug reports and feature requests. Track status of your submissions.",
        "tips": [
            "Use 'Bug' type for issues and 'Feature' for enhancement requests",
            "Provide detailed descriptions to help prioritize",
            "Check 'My Submissions' to see status updates",
            "Admins can triage feedback and update status/priority",
        ],
    },
}


def get_help(page_name: str) -> Optional[Dict[str, Any]]:
    """
    Get help content for a specific page.

    Args:
        page_name: Name of the page

    Returns:
        Help dict with title, description, and tips, or None if not found
    """
    return HELP_TOPICS.get(page_name)


def search_help(query: str) -> List[Dict[str, Any]]:
    """
    Search help topics by fuzzy matching on title and description.

    Args:
        query: Search query string

    Returns:
        List of matching help dicts, ordered by relevance
    """
    if not query or not query.strip():
        return []

    query_lower = query.lower()
    results = []

    for page_name, help_data in HELP_TOPICS.items():
        title = help_data.get("title", "").lower()
        description = help_data.get("description", "").lower()

        # Simple fuzzy matching: check if query appears in title or description
        score = 0
        if query_lower in title:
            score += 10  # Title matches are more important
        if query_lower in description:
            score += 5
        if query_lower in page_name.lower():
            score += 3

        # Check tips
        for tip in help_data.get("tips", []):
            if query_lower in tip.lower():
                score += 1

        if score > 0:
            results.append({
                **help_data,
                "page_name": page_name,
                "score": score,
            })

    # Sort by score descending
    results.sort(key=lambda x: x["score"], reverse=True)
    return results
