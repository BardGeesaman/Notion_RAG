"""Integrate design extraction into repository import workflow."""
from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Experiment
from amprenta_rag.ingestion.design_extraction import (
    detect_design_type,
    extract_geo_design,
    extract_mw_design,
    extract_sample_groups,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def apply_design_extraction(
    experiment_id: UUID,
    repository: str,
    sample_names: Optional[List[str]] = None,
    sample_attributes: Optional[Dict[str, List[str]]] = None,
    study_description: Optional[str] = None,
    factors: Optional[Dict[str, Any]] = None,
) -> bool:
    """
    Apply design extraction to an experiment and update its design fields.

    Args:
        experiment_id: UUID of the experiment to update
        repository: Source repository (GEO, MW, etc.)
        sample_names: List of sample names
        sample_attributes: Sample attribute mappings
        study_description: Study description text
        factors: MW-style factor definitions

    Returns:
        True if design was extracted and saved, False otherwise
    """
    db = next(get_db())
    try:
        experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
        if not experiment:
            logger.warning("[DESIGN] Experiment %s not found", experiment_id)
            return False

        if repository.upper() in ("GEO",):
            result = extract_geo_design({
                "Samples": sample_names or [],
                "Attributes": sample_attributes or {},
                "Summary": study_description or "",
            })
        elif repository.upper() in ("MW", "MW_METABOLOMICS", "MW_LIPIDOMICS"):
            result = extract_mw_design(factors or {})
        else:
            design_type, confidence = detect_design_type(
                sample_names or [],
                sample_attributes,
                study_description,
            )
            groups = extract_sample_groups(sample_names or [], sample_attributes, design_type)
            result = {
                "design_type": design_type,
                "confidence": confidence,
                "sample_groups": groups,
                "design_metadata": {},
            }

        experiment.design_type = result.get("design_type")
        experiment.sample_groups = result.get("sample_groups")
        experiment.design_metadata = {
            "extraction_confidence": result.get("confidence"),
            "source_repository": repository,
            **(result.get("design_metadata") or {}),
        }

        db.commit()
        logger.info(
            "[DESIGN] Applied design_type=%s (confidence=%.2f) to experiment %s",
            experiment.design_type,
            result.get("confidence", 0),
            experiment_id,
        )
        return True

    except Exception as e:
        logger.error("[DESIGN] Error applying design extraction: %r", e)
        db.rollback()
        return False
    finally:
        db.close()


__all__ = ["apply_design_extraction"]


def batch_apply_design_extraction(repository_filter: Optional[str] = None) -> Dict[str, int]:
    """Apply design extraction to all experiments missing design_type."""
    db = next(get_db())
    try:
        query = db.query(Experiment).filter(Experiment.design_type.is_(None))
        if repository_filter:
            # Filter by repository in design_metadata if available
            query = query.filter(Experiment.design_metadata["source_repository"].astext == repository_filter)

        experiments = query.all()
        results = {"processed": 0, "updated": 0, "failed": 0}

        for exp in experiments:
            results["processed"] += 1
            design_type, confidence = detect_design_type(
                sample_names=[],
                sample_attributes=None,
                study_description=exp.description,
            )
            if confidence > 0.5:
                exp.design_type = design_type
                meta = exp.design_metadata or {}
                meta.update({"confidence": confidence, "auto_detected": True})
                exp.design_metadata = meta
                results["updated"] += 1

        db.commit()
        return results
    except Exception as e:
        db.rollback()
        return {"error": str(e)}
    finally:
        db.close()
