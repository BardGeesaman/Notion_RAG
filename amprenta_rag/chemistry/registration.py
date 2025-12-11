"""Compound registration and duplicate checking utilities."""
from __future__ import annotations

import sqlite3
from pathlib import Path
from datetime import datetime
from typing import Optional

from amprenta_rag.chemistry.database import get_chemistry_db_path
from amprenta_rag.chemistry.normalization import normalize_smiles, compute_molecular_descriptors
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _get_db_path(db_path: Optional[Path] = None) -> Path:
    return db_path or get_chemistry_db_path()


def generate_corporate_id(db_path: Optional[Path] = None) -> str:
    """Generate the next corporate ID (e.g., AMP-00001) using the compound_sequence table."""
    path = _get_db_path(db_path)
    conn = sqlite3.connect(str(path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS compound_sequence (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                prefix TEXT DEFAULT 'AMP',
                next_number INTEGER DEFAULT 1
            );
            """
        )
        conn.commit()

        cursor.execute("SELECT id, prefix, next_number FROM compound_sequence ORDER BY id LIMIT 1;")
        row = cursor.fetchone()
        if not row:
            cursor.execute("INSERT INTO compound_sequence (prefix, next_number) VALUES (?, ?)", ("AMP", 1))
            conn.commit()
            cursor.execute("SELECT id, prefix, next_number FROM compound_sequence ORDER BY id LIMIT 1;")
            row = cursor.fetchone()

        seq_id, prefix, next_number = row
        corporate_id = f"{prefix}-{next_number:05d}"
        cursor.execute(
            "UPDATE compound_sequence SET next_number = next_number + 1 WHERE id = ?",
            (seq_id,),
        )
        conn.commit()
        return corporate_id
    except Exception as exc:
        conn.rollback()
        logger.error("[CHEMISTRY][REG] Failed to generate corporate ID: %r", exc)
        raise
    finally:
        conn.close()


def check_duplicate(smiles: str, db_path: Optional[Path] = None) -> Optional[str]:
    """Check for an existing compound by normalized SMILES or InChIKey. Returns compound_id if found."""
    if not smiles or not smiles.strip():
        return None

    canonical_smiles, inchi_key, _ = normalize_smiles(smiles)
    path = _get_db_path(db_path)
    conn = sqlite3.connect(str(path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT compound_id FROM compounds
            WHERE (canonical_smiles = ? OR smiles = ?) OR inchi_key = ?
            LIMIT 1;
            """,
            (canonical_smiles, canonical_smiles, inchi_key),
        )
        row = cursor.fetchone()
        return row[0] if row else None
    finally:
        conn.close()


def register_compound(
    name: str,
    smiles: str,
    salt_form: Optional[str] = None,
    batch_number: Optional[str] = None,
    parent_id: Optional[str] = None,
    registered_by: Optional[str] = None,
    db_path: Optional[Path] = None,
) -> Optional[str]:
    """
    Register a compound if not already present.

    Returns compound_id (corporate_id) or existing duplicate compound_id.
    """
    if not smiles or not smiles.strip():
        return None

    path = _get_db_path(db_path)

    # Duplicate check
    existing = check_duplicate(smiles, path)
    if existing:
        logger.info("[CHEMISTRY][REG] Duplicate compound detected, returning existing id %s", existing)
        return existing

    corporate_id = generate_corporate_id(path)
    canonical_smiles, inchi_key, molecular_formula = normalize_smiles(smiles)
    descriptors = compute_molecular_descriptors(canonical_smiles or smiles)

    now = datetime.utcnow()

    conn = sqlite3.connect(str(path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            INSERT INTO compounds (
                compound_id,
                corporate_id,
                smiles,
                inchi_key,
                canonical_smiles,
                molecular_formula,
                molecular_weight,
                logp,
                hbd_count,
                hba_count,
                rotatable_bonds,
                salt_form,
                batch_number,
                parent_compound_id,
                registration_date,
                registered_by,
                duplicate_of_id,
                created_at,
                updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                corporate_id,
                corporate_id,
                canonical_smiles or smiles,
                inchi_key,
                canonical_smiles,
                molecular_formula,
                descriptors.get("molecular_weight"),
                descriptors.get("logp"),
                descriptors.get("hbd_count"),
                descriptors.get("hba_count"),
                descriptors.get("rotatable_bonds"),
                salt_form,
                batch_number,
                parent_id,
                now,
                registered_by,
                None,
                now,
                now,
            ),
        )
        conn.commit()
        logger.info("[CHEMISTRY][REG] Registered new compound %s", corporate_id)
        return corporate_id
    except Exception as exc:
        conn.rollback()
        logger.error("[CHEMISTRY][REG] Failed to register compound: %r", exc)
        raise
    finally:
        conn.close()
