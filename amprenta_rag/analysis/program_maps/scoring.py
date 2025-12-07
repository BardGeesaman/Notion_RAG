"""
Signature scoring functions for program-signature mapping.

Functions for computing signature scores against program datasets.
"""

from collections import defaultdict
from typing import Dict, List, Optional, Set

from amprenta_rag.analysis.program_maps.models import ProgramSignatureScore
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_program_datasets(program_page_id: str) -> List[str]:
    """
    Get all dataset IDs linked to a program (DEPRECATED - Notion-based).

    ⚠️ DEPRECATED: This function uses Notion API calls.
    Use Postgres relationships instead: `program.datasets` from Postgres Program model.

    Searches for datasets via direct relations or through experiments.

    Args:
        program_page_id: Notion page ID of program (with dashes)

    Returns:
        List of dataset page IDs

    Example:
        >>> datasets = get_program_datasets("abc-123-def")
        >>> len(datasets) >= 0
        True
    """
    import warnings

    warnings.warn(
        "get_program_datasets() is deprecated (Notion-based). "
        "Use Postgres relationships: program.datasets from Postgres Program model.",
        DeprecationWarning,
        stacklevel=2,
    )
    from amprenta_rag.query.cross_omics.helpers import extract_relation_ids, fetch_notion_page

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

    Scores each signature against all datasets in the program and aggregates
    the results to compute overall scores, coverage, and omics-specific scores.

    Args:
        program_page_id: Notion page ID of program (with dashes)
        signature_ids: Optional list of signature IDs to score (if None, scores all)
        use_cache: Whether to use feature cache

    Returns:
        List of ProgramSignatureScore objects, sorted by overall score (descending)

    Example:
        >>> scores = compute_program_signature_scores("abc-123-def")
        >>> len(scores) >= 0
        True
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

    # Load signatures (DEPRECATED - Notion sync disabled)
    raise RuntimeError(
        "This function requires Notion sync which is disabled. " "Use Postgres-based signature loading instead."
    )

    # Dead code below (never executed)
    all_signature_pages = []
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
        logger.warning("[ANALYSIS][PROGRAM-MAPS] No signatures found to score")
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
