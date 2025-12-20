"""
Experimental design extraction from repository metadata.

Detects study design type (case_control, time_course, intervention, etc.)
and extracts sample group assignments from GEO and other repository metadata.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Design type patterns
CASE_CONTROL_PATTERNS = [
    r"(?i)\b(case|patient|disease[d]?|affected|tumor|cancer)\b",
    r"(?i)\b(control|healthy|normal|unaffected|benign)\b",
]

TIME_PATTERNS = [
    r"(?i)\b(day|week|hour|month|time)\s*[\d]+",
    r"(?i)\b[Tt][\d]+\b",  # T0, T1, T2
    r"(?i)\b(baseline|pre|post|before|after)\b",
    r"(?i)\b[\d]+\s*(h|hr|hrs|d|days?|wk|weeks?)\b",
]

INTERVENTION_PATTERNS = [
    r"(?i)\b(treated|untreated|treatment|placebo|vehicle)\b",
    r"(?i)\b(drug|compound|inhibitor|agonist|antagonist)\b",
    r"(?i)\b(knockdown|knockout|overexpression|siRNA|shRNA)\b",
]

DOSE_PATTERNS = [
    r"(?i)\b[\d.]+\s*(mg|ug|ng|pg|[munp]M|mM|uM|nM)\b",
    r"(?i)\b(low|medium|high)\s*dose\b",
    r"(?i)\bdose[\s-]*(response|dependent)\b",
]


def detect_design_type(
    sample_names: List[str],
    sample_attributes: Optional[Dict[str, List[str]]] = None,
    study_description: Optional[str] = None,
) -> Tuple[str, float]:
    """
    Detect experimental design type from sample metadata.

    Args:
        sample_names: List of sample names/IDs
        sample_attributes: Dict mapping attribute name to list of values
        study_description: Study abstract or description text

    Returns:
        Tuple of (design_type, confidence_score)
        design_type is one of: case_control, time_course, intervention,
                               dose_response, multi_factorial, observational
    """
    all_text = " ".join(sample_names)
    if sample_attributes:
        for values in sample_attributes.values():
            all_text += " " + " ".join(str(v) for v in values)
    if study_description:
        all_text += " " + study_description

    scores = {
        "case_control": 0.0,
        "time_course": 0.0,
        "intervention": 0.0,
        "dose_response": 0.0,
    }

    # Check for case/control patterns
    case_matches = sum(1 for p in CASE_CONTROL_PATTERNS[:1] if re.search(p, all_text))
    control_matches = sum(1 for p in CASE_CONTROL_PATTERNS[1:] if re.search(p, all_text))
    if case_matches > 0 and control_matches > 0:
        scores["case_control"] = 0.8

    # Check for time patterns
    time_matches = sum(1 for p in TIME_PATTERNS if re.search(p, all_text))
    if time_matches >= 2:
        scores["time_course"] = 0.7 + min(0.2, time_matches * 0.05)

    # Check for intervention patterns
    intervention_matches = sum(1 for p in INTERVENTION_PATTERNS if re.search(p, all_text))
    if intervention_matches >= 1:
        scores["intervention"] = 0.6 + min(0.3, intervention_matches * 0.1)

    # Check for dose patterns
    dose_matches = sum(1 for p in DOSE_PATTERNS if re.search(p, all_text))
    if dose_matches >= 1:
        scores["dose_response"] = 0.7 + min(0.2, dose_matches * 0.1)

    # Find best match
    best_type = max(scores, key=lambda k: scores.get(k, 0.0))
    best_score = scores[best_type]

    if best_score < 0.5:
        return ("observational", 0.5)

    # Check for multi-factorial (multiple high scores)
    high_scores = sum(1 for s in scores.values() if s >= 0.6)
    if high_scores >= 2:
        return ("multi_factorial", 0.7)

    return (best_type, best_score)


def extract_sample_groups(
    sample_names: List[str],
    sample_attributes: Optional[Dict[str, List[str]]] = None,
    design_type: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Extract sample group assignments based on design type.

    Args:
        sample_names: List of sample names
        sample_attributes: Dict mapping attribute name to list of values
        design_type: Detected or known design type

    Returns:
        Dict mapping group name to list of sample names
        e.g., {"control": ["S1", "S2"], "case": ["S3", "S4"]}
    """
    groups: Dict[str, Any] = {}

    if design_type == "case_control":
        groups = {"control": [], "case": [], "unknown": []}
        for i, name in enumerate(sample_names):
            name_lower = name.lower()
            attr_text = ""
            if sample_attributes:
                for vals in sample_attributes.values():
                    if i < len(vals):
                        attr_text += " " + str(vals[i]).lower()

            combined = name_lower + attr_text
            if any(re.search(p, combined) for p in [r"\bcontrol\b", r"\bhealthy\b", r"\bnormal\b"]):
                groups["control"].append(name)
            elif any(re.search(p, combined) for p in [r"\bcase\b", r"\bpatient\b", r"\bdisease\b", r"\btumor\b"]):
                groups["case"].append(name)
            else:
                groups["unknown"].append(name)

    elif design_type == "time_course":
        groups = {"timepoints": {}}
        for i, name in enumerate(sample_names):
            # Try to extract timepoint from name
            match = re.search(r"[Tt](\d+)|(\d+)\s*(h|hr|d|day|wk)", name)
            if match:
                tp = match.group(0)
                if tp not in groups["timepoints"]:
                    groups["timepoints"][tp] = []
                groups["timepoints"][tp].append(name)
            else:
                if "unknown" not in groups:
                    groups["unknown"] = []
                groups["unknown"].append(name)

    elif design_type == "intervention":
        groups = {"treated": [], "untreated": [], "unknown": []}
        for i, name in enumerate(sample_names):
            name_lower = name.lower()
            if any(w in name_lower for w in ["treated", "treatment", "drug"]):
                groups["treated"].append(name)
            elif any(w in name_lower for w in ["untreated", "control", "vehicle", "placebo"]):
                groups["untreated"].append(name)
            else:
                groups["unknown"].append(name)

    else:
        # Default: no grouping
        groups = {"all": sample_names}

    # Remove empty groups
    return {k: v for k, v in groups.items() if v}


def extract_geo_design(
    geo_metadata: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract experimental design from GEO series metadata.

    Args:
        geo_metadata: GEO series metadata dict (from Entrez)

    Returns:
        Dict with design_type, confidence, sample_groups, and design_metadata
    """
    result = {
        "design_type": "observational",
        "confidence": 0.5,
        "sample_groups": {},
        "design_metadata": {},
    }

    # Extract sample names
    sample_names: List[str] = []
    sample_attrs: Dict[str, List[str]] = {}

    # GEO metadata structure varies; try common paths
    samples = geo_metadata.get("Samples", [])
    if isinstance(samples, list):
        for s in samples:
            if isinstance(s, dict):
                name = s.get("Title", s.get("Accession", ""))
                if name:
                    sample_names.append(name)
                # Extract characteristics
                chars = s.get("Characteristics", [])
                if isinstance(chars, list):
                    for char in chars:
                        if ":" in str(char):
                            key, val = str(char).split(":", 1)
                            key = key.strip()
                            if key not in sample_attrs:
                                sample_attrs[key] = []
                            sample_attrs[key].append(val.strip())

    # Get study description
    description = geo_metadata.get("Summary", "") or geo_metadata.get("Title", "")

    if sample_names:
        design_type, confidence = detect_design_type(sample_names, sample_attrs, description)
        result["design_type"] = design_type
        result["confidence"] = confidence
        result["sample_groups"] = extract_sample_groups(sample_names, sample_attrs, design_type)
        result["design_metadata"] = {
            "sample_count": len(sample_names),
            "attributes_found": list(sample_attrs.keys()),
        }

    logger.info(
        "[DESIGN] Extracted design_type=%s (confidence=%.2f) from %d samples",
        result["design_type"],
        result["confidence"],
        len(sample_names),
    )

    return result


def extract_mw_design(mw_metadata: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract experimental design from Metabolomics Workbench study metadata.

    MW studies have structured FACTORS sections that define experimental variables.

    Args:
        mw_metadata: MW study metadata dict (from mwTab or API)

    Returns:
        Dict with design_type, confidence, sample_groups, and design_metadata
    """
    result = {
        "design_type": "observational",
        "confidence": 0.5,
        "sample_groups": {},
        "design_metadata": {},
    }

    # Extract factors from MW metadata
    factors = mw_metadata.get("FACTORS", mw_metadata.get("factors", []))
    if isinstance(factors, dict):
        factors = [factors]

    factor_names: List[str] = []
    factor_levels: Dict[str, Any] = {}

    for factor in factors:
        if isinstance(factor, dict):
            name = factor.get("name", factor.get("FACTOR_NAME", ""))
            levels = factor.get("levels", factor.get("FACTOR_LEVELS", []))
            if name:
                factor_names.append(name.lower())
                factor_levels[name] = levels if isinstance(levels, list) else [levels]

    # Detect design type from factor names
    factor_text = " ".join(factor_names)

    if any(w in factor_text for w in ["disease", "diagnosis", "condition", "status", "case", "control"]):
        result["design_type"] = "case_control"
        result["confidence"] = 0.85
    elif any(w in factor_text for w in ["time", "timepoint", "day", "hour", "week"]):
        result["design_type"] = "time_course"
        result["confidence"] = 0.85
    elif any(w in factor_text for w in ["treatment", "drug", "compound", "intervention"]):
        result["design_type"] = "intervention"
        result["confidence"] = 0.80
    elif any(w in factor_text for w in ["dose", "concentration"]):
        result["design_type"] = "dose_response"
        result["confidence"] = 0.85
    elif len(factor_names) >= 2:
        result["design_type"] = "multi_factorial"
        result["confidence"] = 0.75

    # Extract sample groups from subject data
    subjects = mw_metadata.get("SUBJECTS", mw_metadata.get("subjects", []))
    if isinstance(subjects, list) and subjects:
        groups: Dict[str, List[str]] = {}
        for subj in subjects:
            if isinstance(subj, dict):
                subj_id = subj.get("SUBJECT_ID", subj.get("subject_id", ""))
                # Try to get factor value for grouping
                for fname in factor_names[:1]:  # Use first factor for primary grouping
                    fval = subj.get(fname, subj.get(fname.upper(), "unknown"))
                    fval_str = str(fval).lower().strip()
                    if fval_str not in groups:
                        groups[fval_str] = []
                    if subj_id:
                        groups[fval_str].append(subj_id)
                    break
                else:
                    if "unknown" not in groups:
                        groups["unknown"] = []
                    if subj_id:
                        groups["unknown"].append(subj_id)

        result["sample_groups"] = {k: v for k, v in groups.items() if v}

    result["design_metadata"] = {
        "factors": factor_names,
        "factor_levels": factor_levels,
        "subject_count": len(subjects) if isinstance(subjects, list) else 0,
    }

    logger.info(
        "[DESIGN] MW extracted design_type=%s (confidence=%.2f) from factors=%s",
        result["design_type"],
        result["confidence"],
        factor_names,
    )

    return result

