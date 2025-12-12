"""Experimental design recommendation engine."""
from __future__ import annotations

from typing import List, Dict, Any

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def recommend_design(
    research_question: str,
    sample_count: int,
    variables: List[str],
) -> Dict[str, Any]:
    """
    Recommend an experimental design based on research question and constraints.
    
    Args:
        research_question: The research question to address
        sample_count: Available or planned sample count
        variables: List of variables/factors to consider
        
    Returns:
        Dict with design_type, rationale, min_samples, considerations
    """
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    variables_str = ", ".join(variables) if variables else "none specified"
    
    prompt = (
        "Analyze this research question and recommend an appropriate experimental design.\n\n"
        f"Research Question: {research_question}\n\n"
        f"Available Samples: {sample_count}\n"
        f"Variables/Factors: {variables_str}\n\n"
        "Common design types:\n"
        "- case_control: Compare two groups (case vs control)\n"
        "- time_course: Measure changes over time\n"
        "- intervention: Before/after treatment comparison\n"
        "- dose_response: Multiple dose levels\n"
        "- factorial: Multiple factors with all combinations\n"
        "- observational: No intervention, observe natural variation\n\n"
        "Respond in JSON format:\n"
        '{"design_type": "case_control", "rationale": "brief explanation", '
        '"min_samples": 40, "considerations": ["consideration1", "consideration2"]}'
    )
    
    try:
        logger.info("[DESIGN] Recommending design for question: %s", research_question[:50])
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3,
            response_format={"type": "json_object"},
        )
        import json
        result = json.loads(resp.choices[0].message.content.strip())  # type: ignore[union-attr]
        
        # Ensure all required fields
        result.setdefault("design_type", "observational")
        result.setdefault("rationale", "Standard observational design")
        result.setdefault("min_samples", sample_count)
        result.setdefault("considerations", [])
        
        logger.info("[DESIGN] Recommended design: %s", result["design_type"])
        return result
    
    except Exception as e:
        logger.error("[DESIGN] Error recommending design: %r", e)
        return {
            "design_type": "observational",
            "rationale": "Error generating recommendation",
            "min_samples": sample_count,
            "considerations": [f"Error: {str(e)}"],
        }


def get_design_requirements(design_type: str) -> Dict[str, Any]:
    """
    Get hardcoded requirements for a design type.
    
    Args:
        design_type: Type of experimental design
        
    Returns:
        Dict with min_samples, min_groups, requirements, description
    """
    requirements_map = {
        "case_control": {
            "min_samples": 20,
            "min_groups": 2,
            "requirements": [
                "Minimum 20 samples per group",
                "Matched controls recommended",
                "Randomization important",
            ],
            "description": "Compare two groups (case vs control)",
        },
        "time_course": {
            "min_samples": 15,
            "min_groups": 3,
            "requirements": [
                "Minimum 3 timepoints",
                "Baseline measurement required",
                "Repeated measures design",
            ],
            "description": "Measure changes over time",
        },
        "intervention": {
            "min_samples": 20,
            "min_groups": 2,
            "requirements": [
                "Minimum 20 samples per group",
                "Pre-intervention baseline",
                "Post-intervention measurement",
            ],
            "description": "Before/after treatment comparison",
        },
        "dose_response": {
            "min_samples": 30,
            "min_groups": 3,
            "requirements": [
                "Minimum 3 dose levels",
                "10 samples per dose level",
                "Include control (dose=0)",
            ],
            "description": "Multiple dose levels to assess response",
        },
        "factorial": {
            "min_samples": 40,
            "min_groups": 4,
            "requirements": [
                "Minimum 2 factors",
                "10 samples per cell (2^k groups)",
                "All factor combinations",
            ],
            "description": "Multiple factors with all combinations",
        },
        "multi_factorial": {
            "min_samples": 60,
            "min_groups": 6,
            "requirements": [
                "Minimum 3 factors",
                "10 samples per cell",
                "Complex interaction analysis",
            ],
            "description": "Three or more factors",
        },
        "observational": {
            "min_samples": 10,
            "min_groups": 1,
            "requirements": [
                "Minimum 10 samples",
                "No intervention",
                "Natural variation study",
            ],
            "description": "Observe natural variation without intervention",
        },
    }
    
    default = {
        "min_samples": 20,
        "min_groups": 2,
        "requirements": ["Standard experimental design"],
        "description": "General experimental design",
    }
    
    return requirements_map.get(design_type.lower(), default)


def validate_design(
    design_type: str,
    sample_count: int,
    group_count: int,
) -> List[str]:
    """
    Validate experimental design against requirements.
    
    Args:
        design_type: Type of experimental design
        sample_count: Total number of samples
        group_count: Number of groups/conditions
        
    Returns:
        List of warning messages if requirements not met
    """
    warnings = []
    
    try:
        requirements = get_design_requirements(design_type)
        min_samples = requirements["min_samples"]
        min_groups = requirements["min_groups"]
        
        # Check sample count
        if sample_count < min_samples:
            warnings.append(
                f"Insufficient samples: {sample_count} < {min_samples} required for {design_type}"
            )
        
        # Check group count
        if group_count < min_groups:
            warnings.append(
                f"Insufficient groups: {group_count} < {min_groups} required for {design_type}"
            )
        
        # Check samples per group
        if group_count > 0:
            samples_per_group = sample_count / group_count
            if samples_per_group < 5:
                warnings.append(
                    f"Very few samples per group: {samples_per_group:.1f} (recommend at least 5-10)"
                )
            elif samples_per_group < 10 and design_type in ["case_control", "intervention", "dose_response"]:
                warnings.append(
                    f"Low samples per group: {samples_per_group:.1f} (recommend at least 10-20 for {design_type})"
                )
        
        logger.debug("[DESIGN] Validation for %s: %d warnings", design_type, len(warnings))
        return warnings
    
    except Exception as e:
        logger.error("[DESIGN] Error validating design: %r", e)
        return [f"Validation error: {str(e)}"]
