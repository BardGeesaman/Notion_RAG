"""
Pathway-aware cross-omics summaries.

Enhances cross-omics summaries with pathway-level insights and enrichment results.
"""

from __future__ import annotations


from amprenta_rag.analysis.enrichment import (
    enrich_dataset_pathways,
    enrich_signature_pathways,
)
from amprenta_rag.logging_utils import get_logger
from typing import Any, Dict, List

logger = get_logger(__name__)


def generate_pathway_aware_dataset_summary(
    dataset_page_id: str,
    top_pathways: int = 10,
    p_value_threshold: float = 0.05,
) -> str:
    """
    Generate a pathway-aware summary for a dataset.

    Args:
        dataset_page_id: Notion page ID of dataset
        top_pathways: Number of top enriched pathways to include
        p_value_threshold: P-value threshold for significance

    Returns:
        Markdown-formatted summary with pathway insights
    """
    logger.info(
        "[ANALYSIS][PATHWAY-SUMMARY] Generating pathway-aware summary for dataset %s",
        dataset_page_id,
    )

    # Perform pathway enrichment
    enrichment_results = enrich_dataset_pathways(
        dataset_page_id=dataset_page_id,
        p_value_threshold=p_value_threshold,
    )

    if not enrichment_results:
        return "No significantly enriched pathways found for this dataset."

    # Build summary
    summary = "## Pathway Enrichment Analysis\n\n"
    summary += f"Found **{len(enrichment_results)}** significantly enriched pathway(s) "
    summary += f"(p < {p_value_threshold}).\n\n"

    # Top pathways
    top_results = enrichment_results[:top_pathways]
    summary += f"### Top {len(top_results)} Enriched Pathways\n\n"

    for i, result in enumerate(top_results, 1):
        pathway = result.pathway
        summary += f"#### {i}. {pathway.name} ({pathway.source})\n\n"
        summary += f"- **Pathway ID**: `{pathway.pathway_id}`\n"
        summary += f"- **P-value**: {result.p_value:.4f}\n"
        summary += f"- **Adjusted P-value**: {result.adjusted_p_value:.4f}\n"
        summary += f"- **Enrichment Ratio**: {result.enrichment_ratio:.2f}x\n"
        summary += f"- **Matched Features**: {result.input_features}/{result.pathway_size}\n"

        if result.matched_features:
            summary += f"- **Features**: {', '.join(result.matched_features[:10])}"
            if len(result.matched_features) > 10:
                summary += f" (and {len(result.matched_features) - 10} more)"
            summary += "\n"

        if pathway.description:
            summary += f"- **Description**: {pathway.description[:200]}...\n"

        summary += "\n"

    # Cross-omics pathway convergence
    pathway_by_source: Dict[str, List[Any]] = {}
    for result in enrichment_results:
        source = result.pathway.source
        if source not in pathway_by_source:
            pathway_by_source[source] = []
        pathway_by_source[source].append(result)

    if len(pathway_by_source) > 1:
        summary += "### Cross-Database Pathway Convergence\n\n"
        summary += "Pathways found across multiple databases:\n"
        for source, results in pathway_by_source.items():
            summary += f"- **{source}**: {len(results)} pathway(s)\n"
        summary += "\n"

    return summary


def generate_pathway_aware_signature_summary(
    signature_page_id: str,
    top_pathways: int = 10,
    p_value_threshold: float = 0.05,
) -> str:
    """
    Generate a pathway-aware summary for a signature.

    Args:
        signature_page_id: Notion page ID of signature
        top_pathways: Number of top enriched pathways to include
        p_value_threshold: P-value threshold for significance

    Returns:
        Markdown-formatted summary with pathway insights
    """
    logger.info(
        "[ANALYSIS][PATHWAY-SUMMARY] Generating pathway-aware summary for signature %s",
        signature_page_id,
    )

    # Perform pathway enrichment
    enrichment_results = enrich_signature_pathways(
        signature_page_id=signature_page_id,
        p_value_threshold=p_value_threshold,
    )

    if not enrichment_results:
        return "No significantly enriched pathways found for this signature."

    # Build summary
    summary = "## Pathway Enrichment Analysis\n\n"
    summary += f"Found **{len(enrichment_results)}** significantly enriched pathway(s) "
    summary += f"(p < {p_value_threshold}).\n\n"

    # Top pathways
    top_results = enrichment_results[:top_pathways]
    summary += f"### Top {len(top_results)} Enriched Pathways\n\n"

    for i, result in enumerate(top_results, 1):
        pathway = result.pathway
        summary += f"#### {i}. {pathway.name} ({pathway.source})\n\n"
        summary += f"- **Pathway ID**: `{pathway.pathway_id}`\n"
        summary += f"- **P-value**: {result.p_value:.4f}\n"
        summary += f"- **Adjusted P-value**: {result.adjusted_p_value:.4f}\n"
        summary += f"- **Enrichment Ratio**: {result.enrichment_ratio:.2f}x\n"
        summary += f"- **Matched Features**: {result.input_features}/{result.pathway_size}\n"

        if result.matched_features:
            summary += f"- **Features**: {', '.join(result.matched_features[:10])}"
            if len(result.matched_features) > 10:
                summary += f" (and {len(result.matched_features) - 10} more)"
            summary += "\n"

        if pathway.description:
            summary += f"- **Description**: {pathway.description[:200]}...\n"

        summary += "\n"

    # Multi-omics pathway coverage
    feature_types = set()
    for result in enrichment_results:
        feature_types.update(result.pathway.feature_types)

    if len(feature_types) > 1:
        summary += "### Multi-Omics Pathway Coverage\n\n"
        summary += "This signature spans pathways across multiple omics types: "
        summary += f"{', '.join(sorted(feature_types))}\n\n"

    return summary


def integrate_pathway_summary_into_cross_omics(
    base_summary: str,
    pathway_summary: str,
) -> str:
    """
    Integrate pathway summary into a cross-omics summary.

    Args:
        base_summary: Base cross-omics summary
        pathway_summary: Pathway enrichment summary

    Returns:
        Combined summary with pathway insights
    """
    if not pathway_summary or pathway_summary == "No significantly enriched pathways found.":
        return base_summary

    # Append pathway summary to base summary
    combined = base_summary
    if not combined.endswith("\n\n"):
        combined += "\n\n"
    combined += pathway_summary

    return combined

