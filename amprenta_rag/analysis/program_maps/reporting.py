"""
Reporting and Notion integration for program-signature mapping.

Functions for generating reports and updating Notion with program maps.
"""

from typing import Any, Dict

from amprenta_rag.analysis.program_maps.convergence import compute_convergence_indicators
from amprenta_rag.analysis.program_maps.coverage import compute_program_omics_coverage
from amprenta_rag.analysis.program_maps.models import ProgramSignatureMap
from amprenta_rag.analysis.program_maps.scoring import compute_program_signature_scores
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def generate_program_signature_map(
    program_page_id: str,
    top_n: int = 10,
    use_cache: bool = True,
) -> ProgramSignatureMap:
    """
    Generate a complete program-signature mapping.

    Orchestrates the computation of signature scores, omics coverage,
    and convergence indicators to create a complete program map.

    Args:
        program_page_id: Notion page ID of program (with dashes)
        top_n: Number of top signatures to include
        use_cache: Whether to use feature cache

    Returns:
        ProgramSignatureMap object

    Example:
        >>> program_map = generate_program_signature_map("abc-123-def")
        >>> program_map.program_id == "abc-123-def"
        True
    """
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Generating signature map for program %s",
        program_page_id,
    )

    # Compute signature scores
    signature_scores = compute_program_signature_scores(
        program_page_id,
        use_cache=use_cache,
    )

    # Compute omics coverage
    omics_coverage = compute_program_omics_coverage(
        program_page_id,
        use_cache=use_cache,
    )

    # Get top signatures
    top_signatures = signature_scores[:top_n] if signature_scores else []

    # Compute convergence indicators
    convergence_indicators = compute_convergence_indicators(signature_scores)

    return ProgramSignatureMap(
        program_id=program_page_id,
        program_name=omics_coverage.program_name,
        signature_scores=signature_scores,
        omics_coverage=omics_coverage,
        top_signatures=top_signatures,
        convergence_indicators=convergence_indicators,
    )


def generate_program_map_report(
    program_map: ProgramSignatureMap,
    include_all_signatures: bool = False,
) -> str:
    """
    Generate a human-readable program-signature map report.

    Creates a markdown-formatted report with omics coverage, convergence
    indicators, and signature scores.

    Args:
        program_map: ProgramSignatureMap object
        include_all_signatures: Whether to include all signatures or just top N

    Returns:
        Markdown-formatted report string

    Example:
        >>> report = generate_program_map_report(program_map)
        >>> "# Program-Signature Map" in report
        True
    """
    report = f"# Program-Signature Map: {program_map.program_name}\n\n"

    # Omics Coverage
    if program_map.omics_coverage:
        report += "## Omics Coverage\n\n"
        report += f"**Total Datasets**: {program_map.omics_coverage.total_datasets}\n\n"

        if program_map.omics_coverage.datasets_by_omics:
            report += "### Datasets by Omics Type\n\n"
            for omics_type, count in sorted(program_map.omics_coverage.datasets_by_omics.items()):
                features_count = program_map.omics_coverage.features_by_omics.get(omics_type, 0)
                report += f"- **{omics_type.capitalize()}**: {count} dataset(s), {features_count} unique features\n"
            report += "\n"

        report += f"**Coverage Summary**: {program_map.omics_coverage.coverage_summary}\n\n"

    # Convergence Indicators
    if program_map.convergence_indicators:
        report += "## Cross-Omics Convergence\n\n"
        report += (
            f"**Multi-Omics Signatures**: {program_map.convergence_indicators.get('multi_omics_signature_count', 0)}\n"
        )
        report += (
            f"**Convergence Fraction**: {program_map.convergence_indicators.get('convergence_fraction', 0.0):.3f}\n"
        )
        report += f"**Avg Omics per Signature**: {program_map.convergence_indicators.get('avg_omics_per_signature', 0.0):.2f}\n\n"

    # Top Signatures
    signatures_to_show = program_map.signature_scores if include_all_signatures else program_map.top_signatures

    report += "## Signature Scores\n\n"
    report += f"**Total Matching Signatures**: {len(program_map.signature_scores)}\n\n"

    if signatures_to_show:
        report += f"### Top {len(signatures_to_show)} Signatures\n\n"
        for i, score in enumerate(signatures_to_show, 1):
            report += f"#### {i}. {score.signature_name}\n\n"
            report += f"- **Overall Score**: {score.overall_score:.3f}\n"
            report += f"- **Coverage**: {score.coverage_fraction:.1%} ({len(score.matching_datasets)}/{program_map.omics_coverage.total_datasets if program_map.omics_coverage else 0} datasets)\n"

            if score.score_by_omics:
                report += "- **Scores by Omics**:\n"
                for omics_type, sig_score in sorted(score.score_by_omics.items()):
                    report += f"  - {omics_type.capitalize()}: {sig_score:.3f}\n"

            report += "\n"
    else:
        report += "No matching signatures found.\n\n"

    return report


def update_notion_with_program_map(program_map: ProgramSignatureMap) -> None:
    """
    Update Notion program page with program-signature map summary (DEPRECATED).

    ⚠️ DEPRECATED: This function uses Notion API calls.
    Consider migrating to Postgres-based program updates instead.

    Attempts to update the program page with omics coverage and top signatures
    information. Searches for suitable properties and updates them if found.

    Args:
        program_map: ProgramSignatureMap object

    Example:
        >>> update_notion_with_program_map(program_map)
    """
    import warnings

    warnings.warn(
        "update_notion_with_program_map() is deprecated. " "Consider using Postgres-based program updates instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    import requests

    from amprenta_rag.clients.notion_client import notion_headers
    from amprenta_rag.config import get_config

    cfg = get_config()

    try:
        # Fetch current page to check for summary properties
        url = f"{cfg.notion.base_url}/pages/{program_map.program_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        page = resp.json()
        props = page.get("properties", {}) or {}

        # Try to find suitable properties for summary
        updates: Dict[str, Any] = {}

        # Update omics coverage if property exists
        if program_map.omics_coverage:
            coverage_prop_candidates = ["Omics Coverage", "Coverage Summary", "Summary"]
            for candidate in coverage_prop_candidates:
                if candidate in props:
                    prop_data = props[candidate]
                    if prop_data.get("type") in ["rich_text", "text"]:
                        coverage_text = program_map.omics_coverage.coverage_summary
                        updates[candidate] = {
                            "rich_text": [
                                {
                                    "type": "text",
                                    "text": {"content": coverage_text[:2000]},  # Notion limit
                                }
                            ],
                        }
                        break

        # Update top signatures if property exists
        if program_map.top_signatures:
            top_sigs_prop_candidates = ["Top Signatures", "Signature Map Summary", "Summary"]
            for candidate in top_sigs_prop_candidates:
                if candidate in props:
                    prop_data = props[candidate]
                    if prop_data.get("type") in ["rich_text", "text"]:
                        top_sigs_text = "\n".join(
                            f"{i+1}. {sig.signature_name} (score: {sig.overall_score:.3f}, coverage: {sig.coverage_fraction:.1%})"
                            for i, sig in enumerate(program_map.top_signatures[:5])
                        )
                        updates[candidate] = {
                            "rich_text": [
                                {
                                    "type": "text",
                                    "text": {"content": top_sigs_text[:2000]},
                                }
                            ],
                        }
                        break

        if updates:
            payload = {"properties": updates}
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            logger.info(
                "[ANALYSIS][PROGRAM-MAPS] Updated Notion program page %s with map summary",
                program_map.program_id,
            )
        else:
            logger.debug(
                "[ANALYSIS][PROGRAM-MAPS] No suitable properties found for program %s. "
                "Consider adding 'Omics Coverage' or 'Top Signatures' properties.",
                program_map.program_id,
            )

    except Exception as e:
        logger.warning(
            "[ANALYSIS][PROGRAM-MAPS] Error updating Notion with program map for %s: %r",
            program_map.program_id,
            e,
        )
        # Non-blocking - continue even if Notion update fails
