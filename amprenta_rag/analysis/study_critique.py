"""Automated study quality critique service."""
from __future__ import annotations

from typing import List, Dict, Any, Protocol
from uuid import UUID

from amprenta_rag.database.models import Experiment
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def assess_study_quality(experiment_id: UUID, db) -> Dict[str, Any]:
    """
    Assess overall study quality for an experiment.

    Args:
        experiment_id: UUID of the experiment to assess
        db: Database session

    Returns:
        Dictionary with quality_score, issues, and summary
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        return {
            "quality_score": 0,
            "issues": [{"flaw": "Experiment not found", "severity": "high", "suggestion": "Verify experiment ID"}],
            "summary": "Experiment not found",
        }

    # Run all checks
    design_flaws = detect_design_flaws(experiment, db)
    data_gaps = identify_data_gaps(experiment, db)

    # Combine issues
    all_issues = design_flaws.copy()
    for gap in data_gaps:
        all_issues.append({
            "flaw": gap,
            "severity": "medium",
            "suggestion": "Add missing information",
        })

    # Calculate quality score
    high_issues = sum(1 for i in all_issues if i.get("severity") == "high")
    medium_issues = sum(1 for i in all_issues if i.get("severity") == "medium")
    low_issues = sum(1 for i in all_issues if i.get("severity") == "low")

    quality_score = max(0, 100 - (high_issues * 20) - (medium_issues * 10) - (low_issues * 5))

    # Generate summary
    if quality_score >= 80:
        summary = "High quality study with minimal issues"
    elif quality_score >= 60:
        summary = "Moderate quality study with some issues to address"
    elif quality_score >= 40:
        summary = "Low quality study with significant issues"
    else:
        summary = "Very low quality study requiring major improvements"

    return {
        "quality_score": quality_score,
        "issues": all_issues,
        "summary": summary,
    }


def detect_design_flaws(experiment: Any, db) -> List[Dict[str, str]]:
    """
    Detect design flaws in an experiment.

    Args:
        experiment: Experiment object
        db: Database session

    Returns:
        List of flaw dictionaries with flaw, severity, and suggestion
    """
    flaws = []

    # Check sample_groups
    if not experiment.sample_groups or experiment.sample_groups == {}:
        flaws.append({
            "flaw": "Missing sample groups",
            "severity": "high",
            "suggestion": "Define sample groups with group names and sample IDs",
        })
    else:
        # Check sample sizes
        for group_name, group_data in experiment.sample_groups.items():
            if isinstance(group_data, dict):
                samples = group_data.get("samples", [])
                if isinstance(samples, list):
                    n = len(samples)
                    if n < 3:
                        flaws.append({
                            "flaw": f"Low sample size in group '{group_name}' (n={n})",
                            "severity": "high" if n == 0 else "medium",
                            "suggestion": f"Increase sample size in group '{group_name}' to at least 3-5 samples",
                        })
            elif isinstance(group_data, list):
                n = len(group_data)
                if n < 3:
                    flaws.append({
                        "flaw": f"Low sample size in group '{group_name}' (n={n})",
                        "severity": "high" if n == 0 else "medium",
                        "suggestion": f"Increase sample size in group '{group_name}' to at least 3-5 samples",
                    })

        # Check for control group
        has_control = False
        control_keywords = ["control", "ctrl", "reference", "baseline", "wildtype", "wt", "vehicle"]
        for group_name in experiment.sample_groups.keys():
            if any(keyword in str(group_name).lower() for keyword in control_keywords):
                has_control = True
                break

        if not has_control:
            flaws.append({
                "flaw": "Missing control group",
                "severity": "high",
                "suggestion": "Add a control group (e.g., 'control', 'vehicle', 'wildtype') for comparison",
            })

    # Check design_type
    if not experiment.design_type:
        flaws.append({
            "flaw": "Design type not specified",
            "severity": "medium",
            "suggestion": "Specify design type (e.g., case_control, time_course, intervention)",
        })

    return flaws


def identify_data_gaps(experiment: Any, db) -> List[str]:
    """
    Identify missing data fields for an experiment.

    Args:
        experiment: Experiment object
        db: Database session

    Returns:
        List of gap description strings
    """
    gaps = []

    # Check organism
    if not experiment.organism or (isinstance(experiment.organism, list) and len(experiment.organism) == 0):
        gaps.append("Missing organism field")

    # Check description
    if not experiment.description or not experiment.description.strip():
        gaps.append("Missing description")

    # Check linked datasets
    if not experiment.datasets or len(experiment.datasets) == 0:
        gaps.append("No linked datasets")

    return gaps
