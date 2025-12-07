"""
Automatic linking helpers for datasets → programs/experiments.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Tuple

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Experiment as ExperimentModel
from amprenta_rag.database.models import Program as ProgramModel
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _normalize_list(values: Optional[Iterable[str]]) -> List[str]:
    return [v.lower() for v in values or [] if v]


def infer_program_from_metadata(
    diseases: Optional[Iterable[str]] = None,
    keywords: Optional[Iterable[str]] = None,
    filename: Optional[str] = None,
    min_confidence: float = 0.8,
) -> Tuple[Optional[str], float]:
    """
    Infer a Program by matching disease overlap and keyword/name hints.

    Returns (program_id, confidence) or (None, 0.0) if ambiguous/low confidence.
    """
    diseases_norm = set(_normalize_list(diseases))
    keywords_norm = set(_normalize_list(keywords))
    if filename:
        for part in filename.replace("_", " ").replace("-", " ").split():
            keywords_norm.add(part.lower())

    db = next(get_db())
    try:
        candidates: List[Tuple[str, float]] = []
        programs = db.query(ProgramModel).all()

        for program in programs:
            score = 0.0
            prog_disease = set(_normalize_list(program.disease))
            if diseases_norm and prog_disease:
                overlap = diseases_norm & prog_disease
                if overlap:
                    score += 0.6 + 0.1 * len(overlap)

            name_tokens = set(_normalize_list(program.name.split())) if program.name else set()
            if keywords_norm and name_tokens and keywords_norm & name_tokens:
                score += 0.3

            if score > 0:
                candidates.append((str(program.id), score))

        if not candidates:
            return None, 0.0

        # Pick top; ensure no ambiguity
        candidates.sort(key=lambda x: x[1], reverse=True)
        top_id, top_score = candidates[0]
        if top_score < min_confidence:
            return None, top_score
        if len(candidates) > 1 and candidates[1][1] == top_score:
            return None, top_score  # ambiguous
        return top_id, top_score
    except Exception as exc:
        logger.warning("[INGEST][AUTO-LINK] Program inference failed: %r", exc)
        return None, 0.0
    finally:
        db.close()


def infer_experiment_from_metadata(
    diseases: Optional[Iterable[str]] = None,
    matrix: Optional[Iterable[str]] = None,
    model_systems: Optional[Iterable[str]] = None,
    min_confidence: float = 0.8,
) -> Tuple[Optional[str], float]:
    """
    Infer an Experiment by matching disease + matrix/model_system overlaps.

    Returns (experiment_id, confidence) or (None, 0.0) if ambiguous/low confidence.
    """
    diseases_norm = set(_normalize_list(diseases))
    matrix_norm = set(_normalize_list(matrix))
    models_norm = set(_normalize_list(model_systems))

    db = next(get_db())
    try:
        candidates: List[Tuple[str, float]] = []
        experiments = db.query(ExperimentModel).all()

        for exp in experiments:
            score = 0.0

            exp_disease = set(_normalize_list(exp.disease))
            exp_matrix = set(_normalize_list(getattr(exp, "matrix", None)))
            exp_models = set(_normalize_list(getattr(exp, "model_systems", None)))

            if diseases_norm and exp_disease:
                overlap = diseases_norm & exp_disease
                if overlap:
                    score += 0.5 + 0.1 * len(overlap)

            if matrix_norm and exp_matrix:
                if matrix_norm & exp_matrix:
                    score += 0.25

            if models_norm and exp_models:
                if models_norm & exp_models:
                    score += 0.25

            if score > 0:
                candidates.append((str(exp.id), score))

        if not candidates:
            return None, 0.0

        candidates.sort(key=lambda x: x[1], reverse=True)
        top_id, top_score = candidates[0]
        if top_score < min_confidence:
            return None, top_score
        if len(candidates) > 1 and candidates[1][1] == top_score:
            return None, top_score  # ambiguous
        return top_id, top_score
    except Exception as exc:
        logger.warning("[INGEST][AUTO-LINK] Experiment inference failed: %r", exc)
        return None, 0.0
    finally:
        db.close()

