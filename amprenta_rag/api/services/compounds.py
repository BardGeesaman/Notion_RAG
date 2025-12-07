from __future__ import annotations

from typing import List, Optional

import sqlite3

from amprenta_rag.chemistry.database import get_chemistry_db_path


def list_compounds() -> List[dict]:
    """List all compounds."""
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT compound_id, smiles, inchi_key, canonical_smiles,
                   molecular_formula, molecular_weight, logp,
                   hbd_count, hba_count, rotatable_bonds
            FROM compounds
            ORDER BY updated_at DESC
            """
        )
        rows = cursor.fetchall()
        return [
            {
                "compound_id": r[0],
                "smiles": r[1],
                "inchi_key": r[2],
                "canonical_smiles": r[3],
                "molecular_formula": r[4],
                "molecular_weight": r[5],
                "logp": r[6],
                "hbd_count": r[7],
                "hba_count": r[8],
                "rotatable_bonds": r[9],
            }
            for r in rows
        ]
    finally:
        conn.close()


def get_compound_by_id(compound_id: str) -> Optional[dict]:
    """Get compound details by ID."""
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT compound_id, smiles, inchi_key, canonical_smiles,
                   molecular_formula, molecular_weight, logp,
                   hbd_count, hba_count, rotatable_bonds
            FROM compounds
            WHERE compound_id = ?
            """,
            (compound_id,),
        )
        r = cursor.fetchone()
        if not r:
            return None
        return {
            "compound_id": r[0],
            "smiles": r[1],
            "inchi_key": r[2],
            "canonical_smiles": r[3],
            "molecular_formula": r[4],
            "molecular_weight": r[5],
            "logp": r[6],
            "hbd_count": r[7],
            "hba_count": r[8],
            "rotatable_bonds": r[9],
        }
    finally:
        conn.close()


def get_compound_programs(compound_id: str) -> List[dict]:
    """Return programs linked to a compound (compound_program table)."""
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT program_id, status, notes, created_at
            FROM compound_program
            WHERE compound_id = ?
            """,
            (compound_id,),
        )
        rows = cursor.fetchall()
        return [
            {
                "program_id": r[0],
                "status": r[1],
                "notes": r[2],
                "created_at": r[3],
            }
            for r in rows
        ]
    finally:
        conn.close()

