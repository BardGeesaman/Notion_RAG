"""
Program-Level Multi-Omics Signature Maps.

Computes Program × Signature scoring matrices, coverage analysis,
and convergence indicators for dashboard intelligence.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.signatures.signature_loader import Signature

logger = get_logger(__name__)


@dataclass
class ProgramSignatureScore:
    """
    Represents a signature score for a program.
    
    Attributes:
        program_id: Notion page ID of program
        signature_id: Notion page ID of signature
        program_name: Name of program
        signature_name: Name of signature
        overall_score: Overall signature score (0-1)
        score_by_omics: Scores by omics type
        matching_datasets: List of dataset IDs that match this signature
        coverage_fraction: Fraction of program datasets that match
    """
    program_id: str
    signature_id: str
    program_name: str
    signature_name: str
    overall_score: float
    score_by_omics: Dict[str, float] = field(default_factory=dict)
    matching_datasets: List[str] = field(default_factory=list)
    coverage_fraction: float = 0.0


@dataclass
class ProgramOmicsCoverage:
    """
    Represents omics coverage for a program.
    
    Attributes:
        program_id: Notion page ID of program
        program_name: Name of program
        total_datasets: Total number of datasets in program
        datasets_by_omics: Number of datasets by omics type
        features_by_omics: Number of unique features by omics type
        coverage_summary: Summary of omics coverage
    """
    program_id: str
    program_name: str
    total_datasets: int
    datasets_by_omics: Dict[str, int] = field(default_factory=dict)
    features_by_omics: Dict[str, int] = field(default_factory=dict)
    coverage_summary: str = ""


@dataclass
class ProgramSignatureMap:
    """
    Represents a complete program-signature mapping.
    
    Attributes:
        program_id: Notion page ID of program
        program_name: Name of program
        signature_scores: List of ProgramSignatureScore objects
        omics_coverage: ProgramOmicsCoverage object
        top_signatures: Top N signatures by score
        convergence_indicators: Cross-omics convergence metrics
    """
    program_id: str
    program_name: str
    signature_scores: List[ProgramSignatureScore] = field(default_factory=list)
    omics_coverage: Optional[ProgramOmicsCoverage] = None
    top_signatures: List[ProgramSignatureScore] = field(default_factory=list)
    convergence_indicators: Dict[str, float] = field(default_factory=dict)


def get_program_datasets(program_page_id: str) -> List[str]:
    """
    Get all dataset IDs linked to a program.
    
    Args:
        program_page_id: Notion page ID of program (with dashes)
        
    Returns:
        List of dataset page IDs
    """
    from amprenta_rag.query.cross_omics.helpers import fetch_notion_page, extract_relation_ids
    
    try:
        page = fetch_notion_page(program_page_id)
        props = page.get("properties", {}) or {}
        
        # Try to find datasets via experiments or directly
        dataset_ids: Set[str] = set()
        
        # Check for direct dataset relations
        dataset_relation_candidates = [
            "Related Datasets",
            "Datasets",
            "Experimental Data Assets",
            "Program Datasets",
        ]
        
        for candidate in dataset_relation_candidates:
            if candidate in props:
                dataset_ids.update(extract_relation_ids(page, candidate))
        
        # Also check via experiments
        experiment_relation_candidates = [
            "Related Experiments",
            "Experiments",
        ]
        
        for candidate in experiment_relation_candidates:
            if candidate in props:
                experiment_ids = extract_relation_ids(page, candidate)
                # Fetch each experiment and get its datasets
                for exp_id in experiment_ids:
                    try:
                        exp_page = fetch_notion_page(exp_id)
                        if exp_page:
                            exp_dataset_candidates = [
                                "Related Datasets",
                                "Datasets",
                            ]
                            for exp_candidate in exp_dataset_candidates:
                                dataset_ids.update(extract_relation_ids(exp_page, exp_candidate))
                    except Exception as e:
                        logger.debug(
                            "[ANALYSIS][PROGRAM-MAPS] Error fetching experiment %s: %r",
                            exp_id,
                            e,
                        )
                        continue
        
        return list(dataset_ids)
        
    except Exception as e:
        logger.error(
            "[ANALYSIS][PROGRAM-MAPS] Error getting datasets for program %s: %r",
            program_page_id,
            e,
        )
        return []


def compute_program_signature_scores(
    program_page_id: str,
    signature_ids: Optional[List[str]] = None,
    use_cache: bool = True,
) -> List[ProgramSignatureScore]:
    """
    Compute signature scores for all datasets in a program.
    
    Args:
        program_page_id: Notion page ID of program (with dashes)
        signature_ids: Optional list of signature IDs to score (if None, scores all)
        use_cache: Whether to use feature cache
        
    Returns:
        List of ProgramSignatureScore objects
    """
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Computing signature scores for program %s",
        program_page_id,
    )
    
    # Get program name
    program_name = f"Program {program_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import fetch_notion_page
        page = fetch_notion_page(program_page_id)
        props = page.get("properties", {}) or {}
        name_prop = props.get("Program", {}).get("title", []) or []
        if not name_prop:
            name_prop = props.get("Name", {}).get("title", []) or []
        if name_prop:
            program_name = name_prop[0].get("plain_text", program_name)
    except Exception as e:
        logger.debug("[ANALYSIS][PROGRAM-MAPS] Could not fetch program name: %r", e)
    
    # Get all datasets in program
    dataset_ids = get_program_datasets(program_page_id)
    
    if not dataset_ids:
        logger.warning(
            "[ANALYSIS][PROGRAM-MAPS] No datasets found for program %s",
            program_page_id,
        )
        return []
    
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Found %d datasets for program %s",
        len(dataset_ids),
        program_page_id,
    )
    
    # Load signatures
    from amprenta_rag.ingestion.signature_matching.signature_loader import (
        fetch_all_signatures_from_notion,
        load_signature_from_notion_page,
    )
    
    all_signature_pages = fetch_all_signatures_from_notion()
    all_signatures = []
    sig_name_to_id: Dict[str, str] = {}
    
    for sig_page in all_signature_pages:
        sig_id = sig_page.get("id", "")
        props = sig_page.get("properties", {}) or {}
        name_prop = props.get("Name", {}).get("title", []) or []
        sig_name = name_prop[0].get("plain_text", "") if name_prop else ""
        
        # Filter by signature_ids if provided
        if signature_ids and sig_id not in signature_ids:
            continue
        
        sig = load_signature_from_notion_page(sig_page)
        if sig:
            all_signatures.append(sig)
            if sig_name:
                sig_name_to_id[sig_name] = sig_id
    
    if not all_signatures:
        logger.warning(
            "[ANALYSIS][PROGRAM-MAPS] No signatures found to score"
        )
        return []
    
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Scoring %d signatures against %d datasets",
        len(all_signatures),
        len(dataset_ids),
    )
    
    # Score each signature against program datasets
    signature_scores: Dict[str, ProgramSignatureScore] = {}
    
    for signature in all_signatures:
        matching_datasets: List[str] = []
        scores_by_dataset: Dict[str, float] = {}
        scores_by_omics: Dict[str, List[float]] = defaultdict(list)
        
        for dataset_id in dataset_ids:
            try:
                # Extract features from dataset
                features_by_type = extract_dataset_features_by_type(
                    dataset_id,
                    use_cache=use_cache,
                )
                
                # Score signature against dataset
                from amprenta_rag.ingestion.multi_omics_scoring import score_multi_omics_signature_against_dataset
                
                score_result = score_multi_omics_signature_against_dataset(
                    signature=signature,
                    dataset_features_by_type=features_by_type,
                )
                
                if score_result and score_result.total_score > 0:
                    scores_by_dataset[dataset_id] = score_result.total_score
                    matching_datasets.append(dataset_id)
                    
                    # Track scores by omics type
                    for comp_match in score_result.component_matches:
                        if comp_match.matched:
                            omics_type = comp_match.feature_type
                            scores_by_omics[omics_type].append(comp_match.match_score)
                
            except Exception as e:
                logger.debug(
                    "[ANALYSIS][PROGRAM-MAPS] Error scoring signature %s against dataset %s: %r",
                    signature.name,
                    dataset_id,
                    e,
                )
                continue
        
        if matching_datasets:
            # Compute overall score (average across matching datasets)
            overall_score = sum(scores_by_dataset.values()) / len(scores_by_dataset) if scores_by_dataset else 0.0
            
            # Compute average score by omics type
            score_by_omics: Dict[str, float] = {}
            for omics_type, scores in scores_by_omics.items():
                score_by_omics[omics_type] = sum(scores) / len(scores) if scores else 0.0
            
            # Coverage fraction
            coverage_fraction = len(matching_datasets) / len(dataset_ids) if dataset_ids else 0.0
            
            # Get signature page ID from mapping
            sig_page_id = sig_name_to_id.get(signature.name, signature.name)
            signature_scores[sig_page_id] = ProgramSignatureScore(
                program_id=program_page_id,
                signature_id=sig_page_id,
                program_name=program_name,
                signature_name=signature.name,
                overall_score=overall_score,
                score_by_omics=score_by_omics,
                matching_datasets=matching_datasets,
                coverage_fraction=coverage_fraction,
            )
    
    result = list(signature_scores.values())
    result.sort(key=lambda s: s.overall_score, reverse=True)
    
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Computed scores for %d signatures for program %s",
        len(result),
        program_page_id,
    )
    
    return result


def compute_program_omics_coverage(
    program_page_id: str,
    use_cache: bool = True,
) -> ProgramOmicsCoverage:
    """
    Compute omics coverage for a program.
    
    Args:
        program_page_id: Notion page ID of program (with dashes)
        use_cache: Whether to use feature cache
        
    Returns:
        ProgramOmicsCoverage object
    """
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Computing omics coverage for program %s",
        program_page_id,
    )
    
    # Get program name
    program_name = f"Program {program_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import fetch_notion_page
        page = fetch_notion_page(program_page_id)
        props = page.get("properties", {}) or {}
        name_prop = props.get("Program", {}).get("title", []) or []
        if not name_prop:
            name_prop = props.get("Name", {}).get("title", []) or []
        if name_prop:
            program_name = name_prop[0].get("plain_text", program_name)
    except Exception as e:
        logger.debug("[ANALYSIS][PROGRAM-MAPS] Could not fetch program name: %r", e)
    
    # Get all datasets
    dataset_ids = get_program_datasets(program_page_id)
    
    if not dataset_ids:
        return ProgramOmicsCoverage(
            program_id=program_page_id,
            program_name=program_name,
            total_datasets=0,
        )
    
    # Analyze datasets by omics type
    datasets_by_omics: Dict[str, int] = defaultdict(int)
    all_features_by_omics: Dict[str, Set[str]] = defaultdict(set)
    
    for dataset_id in dataset_ids:
        try:
            features_by_type = extract_dataset_features_by_type(
                dataset_id,
                use_cache=use_cache,
            )
            
            for omics_type, features in features_by_type.items():
                if features:
                    datasets_by_omics[omics_type] += 1
                    all_features_by_omics[omics_type].update(features)
        
        except Exception as e:
            logger.debug(
                "[ANALYSIS][PROGRAM-MAPS] Error analyzing dataset %s: %r",
                dataset_id,
                e,
            )
            continue
    
    features_by_omics = {k: len(v) for k, v in all_features_by_omics.items()}
    
    # Generate coverage summary
    coverage_parts = []
    for omics_type in sorted(datasets_by_omics.keys()):
        datasets_count = datasets_by_omics[omics_type]
        features_count = features_by_omics.get(omics_type, 0)
        coverage_parts.append(
            f"{omics_type.capitalize()}: {datasets_count} dataset(s), {features_count} unique features"
        )
    
    coverage_summary = "; ".join(coverage_parts) if coverage_parts else "No omics data"
    
    return ProgramOmicsCoverage(
        program_id=program_page_id,
        program_name=program_name,
        total_datasets=len(dataset_ids),
        datasets_by_omics=dict(datasets_by_omics),
        features_by_omics=features_by_omics,
        coverage_summary=coverage_summary,
    )


def compute_convergence_indicators(
    signature_scores: List[ProgramSignatureScore],
) -> Dict[str, float]:
    """
    Compute cross-omics convergence indicators.
    
    Args:
        signature_scores: List of ProgramSignatureScore objects
        
    Returns:
        Dictionary of convergence metrics
    """
    if not signature_scores:
        return {}
    
    # Count signatures with multi-omics support
    multi_omics_signatures = 0
    total_omics_types = 0
    
    for score in signature_scores:
        omics_count = len([s for s in score.score_by_omics.values() if s > 0])
        if omics_count > 1:
            multi_omics_signatures += 1
        total_omics_types += omics_count
    
    convergence_fraction = multi_omics_signatures / len(signature_scores) if signature_scores else 0.0
    avg_omics_per_signature = total_omics_types / len(signature_scores) if signature_scores else 0.0
    
    return {
        "convergence_fraction": convergence_fraction,
        "avg_omics_per_signature": avg_omics_per_signature,
        "multi_omics_signature_count": multi_omics_signatures,
    }


def generate_program_signature_map(
    program_page_id: str,
    top_n: int = 10,
    use_cache: bool = True,
) -> ProgramSignatureMap:
    """
    Generate a complete program-signature mapping.
    
    Args:
        program_page_id: Notion page ID of program (with dashes)
        top_n: Number of top signatures to include
        use_cache: Whether to use feature cache
        
    Returns:
        ProgramSignatureMap object
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
    
    Args:
        program_map: ProgramSignatureMap object
        include_all_signatures: Whether to include all signatures or just top N
        
    Returns:
        Markdown-formatted report
    """
    report = f"# Program-Signature Map: {program_map.program_name}\n\n"
    
    # Omics Coverage
    if program_map.omics_coverage:
        report += f"## Omics Coverage\n\n"
        report += f"**Total Datasets**: {program_map.omics_coverage.total_datasets}\n\n"
        
        if program_map.omics_coverage.datasets_by_omics:
            report += f"### Datasets by Omics Type\n\n"
            for omics_type, count in sorted(program_map.omics_coverage.datasets_by_omics.items()):
                features_count = program_map.omics_coverage.features_by_omics.get(omics_type, 0)
                report += f"- **{omics_type.capitalize()}**: {count} dataset(s), {features_count} unique features\n"
            report += "\n"
        
        report += f"**Coverage Summary**: {program_map.omics_coverage.coverage_summary}\n\n"
    
    # Convergence Indicators
    if program_map.convergence_indicators:
        report += f"## Cross-Omics Convergence\n\n"
        report += f"**Multi-Omics Signatures**: {program_map.convergence_indicators.get('multi_omics_signature_count', 0)}\n"
        report += f"**Convergence Fraction**: {program_map.convergence_indicators.get('convergence_fraction', 0.0):.3f}\n"
        report += f"**Avg Omics per Signature**: {program_map.convergence_indicators.get('avg_omics_per_signature', 0.0):.2f}\n\n"
    
    # Top Signatures
    signatures_to_show = program_map.signature_scores if include_all_signatures else program_map.top_signatures
    
    report += f"## Signature Scores\n\n"
    report += f"**Total Matching Signatures**: {len(program_map.signature_scores)}\n\n"
    
    if signatures_to_show:
        report += f"### Top {len(signatures_to_show)} Signatures\n\n"
        for i, score in enumerate(signatures_to_show, 1):
            report += f"#### {i}. {score.signature_name}\n\n"
            report += f"- **Overall Score**: {score.overall_score:.3f}\n"
            report += f"- **Coverage**: {score.coverage_fraction:.1%} ({len(score.matching_datasets)}/{program_map.omics_coverage.total_datasets if program_map.omics_coverage else 0} datasets)\n"
            
            if score.score_by_omics:
                report += f"- **Scores by Omics**:\n"
                for omics_type, sig_score in sorted(score.score_by_omics.items()):
                    report += f"  - {omics_type.capitalize()}: {sig_score:.3f}\n"
            
            report += "\n"
    else:
        report += "No matching signatures found.\n\n"
    
    return report


def update_notion_with_program_map(program_map: ProgramSignatureMap) -> None:
    """
    Update Notion program page with program-signature map summary.
    
    Args:
        program_map: ProgramSignatureMap object
    """
    # DEPRECATED: Notion imports removed
    # from amprenta_rag.clients.notion_client import notion_headers
    from amprenta_rag.config import get_config
    import requests
    
    def notion_headers() -> Dict[str, str]:
        """DEPRECATED: Notion support removed. Returns empty headers dict."""
        logger.debug("[PROGRAM-SIGNATURE-MAPS] notion_headers() deprecated - Notion support removed")
        return {}
    
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
