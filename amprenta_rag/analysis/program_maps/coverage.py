"""
Omics coverage computation for program-signature mapping.

Functions for analyzing omics data coverage across program datasets.
Postgres is now the source of truth - Notion support has been removed.
"""

from collections import defaultdict
from typing import Dict, List, Set
from uuid import UUID

from amprenta_rag.analysis.program_maps.models import ProgramOmicsCoverage
from amprenta_rag.analysis.program_maps.scoring import get_program_datasets
from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Program as ProgramModel
from amprenta_rag.ingestion.multi_omics_scoring import extract_dataset_features_by_type
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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


def compute_program_omics_coverage(
    program_page_id: str,
    use_cache: bool = True,
) -> ProgramOmicsCoverage:
    """
    Compute omics coverage for a program.

    Analyzes all datasets in a program to determine:
    - Number of datasets by omics type
    - Number of unique features by omics type
    - Overall coverage summary

    Args:
        program_page_id: Program ID (UUID or legacy Notion page ID)
        use_cache: Whether to use feature cache

    Returns:
        ProgramOmicsCoverage object

    Example:
        >>> coverage = compute_program_omics_coverage("abc-123-def")
        >>> coverage.total_datasets >= 0
        True
    """
    logger.info(
        "[ANALYSIS][PROGRAM-MAPS] Computing omics coverage for program %s",
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
        return ProgramOmicsCoverage(
            program_id=program_page_id,
            program_name=program_name,
            total_datasets=0,
            datasets_by_omics={},
            features_by_omics={},
            coverage_summary="No datasets found",
        )

    # Analyze datasets by omics type
    datasets_by_omics: Dict[str, int] = defaultdict(int)
    features_by_omics: Dict[str, Set[str]] = defaultdict(set)

    for dataset_id in dataset_ids:
        try:
            features_by_type = extract_dataset_features_by_type(
                dataset_id,
                use_cache=use_cache,
            )

            # Count datasets and collect features by omics type
            for omics_type, features in features_by_type.items():
                if features:
                    datasets_by_omics[omics_type] += 1
                    features_by_omics[omics_type].update(features)

        except Exception as e:
            logger.debug(
                "[ANALYSIS][PROGRAM-MAPS] Error analyzing dataset %s: %r",
                dataset_id,
                e,
            )
            continue

    # Convert feature sets to counts
    features_count_by_omics: Dict[str, int] = {
        omics_type: len(features) for omics_type, features in features_by_omics.items()
    }

    # Generate coverage summary
    coverage_parts: List[str] = []
    for omics_type in sorted(datasets_by_omics.keys()):
        datasets_count = datasets_by_omics[omics_type]
        features_count = features_count_by_omics.get(omics_type, 0)
        coverage_parts.append(
            f"{omics_type.capitalize()}: {datasets_count} dataset(s), {features_count} unique features"
        )

    coverage_summary = "; ".join(coverage_parts) if coverage_parts else "No omics data"

    return ProgramOmicsCoverage(
        program_id=program_page_id,
        program_name=program_name,
        total_datasets=len(dataset_ids),
        datasets_by_omics=dict(datasets_by_omics),
        features_by_omics=features_count_by_omics,
        coverage_summary=coverage_summary,
    )
