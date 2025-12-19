"""
Signature scoring functions for program-signature mapping.

Functions for computing signature scores against program datasets.
Postgres is now the source of truth - Notion support has been removed.
"""

from collections import defaultdict
from typing import Dict, List, Optional
from uuid import UUID

from amprenta_rag.analysis.program_maps.models import ProgramSignatureScore
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Program as ProgramModel
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import load_signatures_from_postgres

logger = get_logger(__name__)


def get_program_datasets(program_page_id: str) -> List[str]:
    """
    Get all dataset IDs linked to a program from Postgres.

    Args:
        program_page_id: Program ID (UUID or legacy Notion page ID)

    Returns:
        List of dataset IDs (as strings)

    Example:
        >>> datasets = get_program_datasets("abc-123-def")
        >>> len(datasets) >= 0
        True
    """
    with db_session() as db:
        try:
            program = None
            try:
                program_uuid = UUID(program_page_id)
                program = db.query(ProgramModel).filter(ProgramModel.id == program_uuid).first()
            except ValueError:
                program = db.query(ProgramModel).filter(ProgramModel.notion_page_id == program_page_id).first()

            if not program:
                return []

            return [str(dataset.id) for dataset in program.datasets]
        except Exception as e:
            logger.warning(
                "[ANALYSIS][PROGRAM-MAPS] Error getting datasets for program %s: %r",
                program_page_id,
                e,
            )
            return []


def _get_program_name(program_page_id: str) -> str:
    """Get program name from Postgres."""
    with db_session() as db:
        try:
            program = None
            try:
                program_uuid = UUID(program_page_id)
                program = db.query(ProgramModel).filter(ProgramModel.id == program_uuid).first()
            except ValueError:
                program = db.query(ProgramModel).filter(ProgramModel.notion_page_id == program_page_id).first()

            if program and program.name:
                return program.name
            return f"Program {program_page_id[:8]}"
        except Exception as e:
            logger.debug("[ANALYSIS][PROGRAM-MAPS] Could not fetch program name: %r", e)
            return f"Program {program_page_id[:8]}"


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
        program_page_id: Program ID (UUID or legacy Notion page ID)
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

    program_name = _get_program_name(program_page_id)

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

    # Load signatures from Postgres
    all_signatures = load_signatures_from_postgres()
    sig_name_to_id: Dict[str, str] = {}

    # Filter by signature_ids if provided
    if signature_ids:
        all_signatures = [s for s in all_signatures if str(getattr(s, 'id', '')) in signature_ids]

    for sig in all_signatures:
        if sig.name:
            sig_name_to_id[sig.name] = str(getattr(sig, 'id', sig.name))

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

            # Get signature ID
            sig_id = sig_name_to_id.get(signature.name, signature.name)
            signature_scores[sig_id] = ProgramSignatureScore(
                program_id=program_page_id,
                signature_id=sig_id,
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
