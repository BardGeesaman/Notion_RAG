from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

from amprenta_rag.chemistry.database import get_chemistry_db_path
from amprenta_rag.chemistry.schema import CompoundSignatureLink
from amprenta_rag.logging_utils import get_logger
import sqlite3

logger = get_logger(__name__)


def _get_conn(db_path: Optional[str] = None) -> sqlite3.Connection:
    path = db_path or get_chemistry_db_path()
    return sqlite3.connect(str(path))


def link_compound_to_program(
    compound_id: str,
    program_id: str,
    role: str = "HIT",
    db_path: Optional[str] = None,
) -> None:
    """
    Link a compound to a program with a role (HIT, LEAD, TOOL, CANDIDATE).
    """
    conn = _get_conn(db_path)
    try:
        conn.execute(
            """
            INSERT OR REPLACE INTO compound_program (compound_id, program_id, role, created_at, updated_at)
            VALUES (?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
            """,
            (compound_id, program_id, role),
        )
        conn.commit()
        logger.info("[CHEMISTRY][LINKING] Linked compound %s to program %s as %s", compound_id, program_id, role)
    except Exception as e:
        conn.rollback()
        logger.error("[CHEMISTRY][LINKING] Error linking compound %s to program %s: %r", compound_id, program_id, e)
        raise
    finally:
        conn.close()


def link_compound_to_signature(
    compound_id: str,
    signature_id: str,
    effect_type: str,
    correlation: Optional[float] = None,
    p_value: Optional[float] = None,
    evidence_source: Optional[str] = None,
    db_path: Optional[str] = None,
) -> None:
    """
    Link a compound to a signature with an effect type.
    effect_type: 'reverses', 'mimics', 'partial', 'unknown'
    """
    conn = _get_conn(db_path)
    try:
        conn.execute(
            """
            INSERT OR REPLACE INTO compound_signature (
                compound_id, signature_id, effect_type, correlation, p_value, evidence_source, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
            """,
            (compound_id, signature_id, effect_type, correlation, p_value, evidence_source),
        )
        conn.commit()
        logger.info(
            "[CHEMISTRY][LINKING] Linked compound %s to signature %s (%s)",
            compound_id,
            signature_id,
            effect_type,
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][LINKING] Error linking compound %s to signature %s: %r",
            compound_id,
            signature_id,
            e,
        )
        raise
    finally:
        conn.close()


def find_compounds_reversing_signature(
    signature_id: str,
    min_correlation: float = -0.5,
    db_path: Optional[str] = None,
) -> List[Tuple[str, Optional[float], Optional[float]]]:
    """
    Return list of (compound_id, correlation, p_value) for compounds with effect_type='reverses'
    filtered by correlation threshold.
    """
    conn = _get_conn(db_path)
    try:
        rows = conn.execute(
            """
            SELECT compound_id, correlation, p_value
            FROM compound_signature
            WHERE signature_id = ? AND effect_type = 'reverses'
            """,
            (signature_id,),
        ).fetchall()
        results = []
        for cid, corr, pval in rows:
            if corr is None or corr <= min_correlation:
                results.append((cid, corr, pval))
        logger.info(
            "[CHEMISTRY][LINKING] Found %d reversing compounds for signature %s",
            len(results),
            signature_id,
        )
        return results
    finally:
        conn.close()


def find_signatures_affected_by_compound(
    compound_id: str,
    db_path: Optional[str] = None,
) -> List[CompoundSignatureLink]:
    """
    Get all signatures linked to a compound.
    """
    conn = _get_conn(db_path)
    try:
        rows = conn.execute(
            """
            SELECT compound_id, signature_id, effect_type, correlation, p_value, evidence_source, created_at, updated_at
            FROM compound_signature
            WHERE compound_id = ?
            """,
            (compound_id,),
        ).fetchall()
        links: List[CompoundSignatureLink] = []
        for row in rows:
            links.append(
                CompoundSignatureLink(
                    compound_id=row[0],
                    signature_id=row[1],
                    effect_type=row[2],
                    correlation=row[3],
                    p_value=row[4],
                    evidence_source=row[5],
                    created_at=row[6],
                    updated_at=row[7],
                )
            )
        logger.info("[CHEMISTRY][LINKING] Found %d signatures linked to compound %s", len(links), compound_id)
        return links
    finally:
        conn.close()


def compute_signature_reversal_score(
    compound_features: Dict[str, float],
    signature_features: Dict[str, float],
) -> Tuple[float, float, str]:
    """
    Given compound effect profile and signature features, compute anti-correlation score.

    Returns (correlation, p_value, effect_type)
    """
    shared = set(compound_features.keys()) & set(signature_features.keys())
    if len(shared) < 3:
        logger.info("[CHEMISTRY][LINKING] Insufficient overlap for reversal score (n=%d)", len(shared))
        return (np.nan, np.nan, "unknown")

    compound_vals = np.array([compound_features[k] for k in shared], dtype=float)
    signature_vals = np.array([signature_features[k] for k in shared], dtype=float)

    # Anti-correlation (reverse) expectation: negative correlation
    corr, pval = stats.pearsonr(compound_vals, signature_vals)
    if np.isnan(corr):
        return (corr, pval, "unknown")

    effect_type = "reverses" if corr < 0 else "mimics"
    return (float(corr), float(pval), effect_type)

