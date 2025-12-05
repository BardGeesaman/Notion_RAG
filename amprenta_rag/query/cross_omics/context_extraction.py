"""
Context extraction for cross-omics reasoning.

Extracts disease, model system, matrix, and other contextual information
from Notion pages to enhance cross-omics summaries.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics.helpers import (
    extract_relation_ids,
    extract_select_values,
    extract_text_property,
    fetch_notion_page,
)

logger = get_logger(__name__)


def extract_disease_context(page: Dict[str, Any]) -> List[str]:
    """
    Extract disease context from a Notion page.
    
    Args:
        page: Notion page dictionary
        
    Returns:
        List of disease names
    """
    diseases = extract_select_values(page, "Disease")
    if not diseases:
        # Try alternative property names
        diseases = extract_select_values(page, "Disease Context")
    return diseases


def extract_matrix_context(page: Dict[str, Any]) -> List[str]:
    """
    Extract matrix context from a Notion page.
    
    Args:
        page: Notion page dictionary
        
    Returns:
        List of matrix types (CSF, plasma, serum, tissue, etc.)
    """
    return extract_select_values(page, "Matrix")


def extract_model_system_context(page: Dict[str, Any]) -> List[str]:
    """
    Extract model system context from a Notion page.
    
    Args:
        page: Notion page dictionary
        
    Returns:
        List of model systems (in vitro, in vivo, patient, cell line, etc.)
    """
    return extract_select_values(page, "Model Systems")


def extract_aggregated_context(
    page_ids: List[str],
    page_type: str = "dataset",
) -> Dict[str, Any]:
    """
    Extract and aggregate context from multiple Notion pages.
    
    Aggregates disease, matrix, and model system information across
    multiple pages to provide comprehensive context.
    
    Args:
        page_ids: List of Notion page IDs
        page_type: Type of pages ("dataset", "experiment", "program")
        
    Returns:
        Dictionary with aggregated context:
        - diseases: Set of all unique diseases
        - matrix: Set of all unique matrix types
        - model_systems: Set of all unique model systems
        - page_count: Number of pages processed
    """
    all_diseases: set[str] = set()
    all_matrix: set[str] = set()
    all_model_systems: set[str] = set()
    processed_count = 0
    
    for page_id in page_ids:
        try:
            page = fetch_notion_page(page_id)
            if not page:
                continue
                
            diseases = extract_disease_context(page)
            matrix = extract_matrix_context(page)
            model_systems = extract_model_system_context(page)
            
            all_diseases.update(diseases)
            all_matrix.update(matrix)
            all_model_systems.update(model_systems)
            processed_count += 1
            
        except Exception as e:
            logger.debug(
                "[RAG][CROSS-OMICS] Error extracting context from %s %s: %r",
                page_type,
                page_id,
                e,
            )
            continue
    
    return {
        "diseases": sorted(list(all_diseases)),
        "matrix": sorted(list(all_matrix)),
        "model_systems": sorted(list(all_model_systems)),
        "page_count": processed_count,
    }


def format_context_for_prompt(context: Dict[str, Any]) -> str:
    """
    Format aggregated context for inclusion in LLM prompts.
    
    Args:
        context: Dictionary from extract_aggregated_context
        
    Returns:
        Formatted string describing the context
    """
    parts: List[str] = []
    
    if context.get("diseases"):
        diseases_str = ", ".join(context["diseases"])
        parts.append(f"Disease context: {diseases_str}")
    
    if context.get("matrix"):
        matrix_str = ", ".join(context["matrix"])
        parts.append(f"Matrix types: {matrix_str}")
    
    if context.get("model_systems"):
        model_str = ", ".join(context["model_systems"])
        parts.append(f"Model systems: {model_str}")
    
    if not parts:
        return "No disease, matrix, or model system context available."
    
    return "\n".join(parts)


def identify_comparative_context(
    context: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    """
    Identify if context suggests comparative analysis (disease vs control).
    
    Looks for patterns that suggest comparative studies:
    - Disease terms present
    - Control/comparison groups mentioned
    - Multiple conditions
    
    Args:
        context: Dictionary from extract_aggregated_context
        
    Returns:
        Dictionary with comparative analysis hints, or None if not applicable
    """
    diseases = context.get("diseases", [])
    if not diseases:
        return None
    
    # Check for control/comparison indicators
    # This is a heuristic - could be enhanced with more sophisticated detection
    comparative_hints: Dict[str, Any] = {
        "has_disease_context": True,
        "diseases": diseases,
        "suggests_comparison": True,
    }
    
    # Check model systems for patient/control patterns
    model_systems = context.get("model_systems", [])
    if any("patient" in ms.lower() or "control" in ms.lower() for ms in model_systems):
        comparative_hints["has_control_group"] = True
    
    return comparative_hints

