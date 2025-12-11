"""Substructure and similarity search utilities using RDKit."""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import List, Optional

from amprenta_rag.chemistry.database import get_chemistry_db_path
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Try RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("[CHEMISTRY][SEARCH] RDKit not available; structure search disabled.")


def _get_db_path(db_path: Optional[Path] = None) -> Path:
    return db_path or get_chemistry_db_path()


def _load_compounds(db_path: Optional[Path] = None) -> List[tuple]:
    path = _get_db_path(db_path)
    conn = sqlite3.connect(str(path))
    try:
        rows = conn.execute(
            "SELECT compound_id, smiles, corporate_id FROM compounds"
        ).fetchall()
        return rows
    finally:
        conn.close()


def substructure_search(query_smarts: str, db_path: Optional[Path] = None) -> List[dict]:
    """Substructure search; returns list of matching compounds."""
    if not RDKIT_AVAILABLE:
        return []
    if not query_smarts:
        return []

    pattern = Chem.MolFromSmarts(query_smarts)
    if pattern is None:
        logger.warning("[CHEMISTRY][SEARCH] Invalid SMARTS: %s", query_smarts)
        return []

    matches = []
    for compound_id, smiles, corporate_id in _load_compounds(db_path):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        if mol.HasSubstructMatch(pattern):
            matches.append(
                {
                    "compound_id": compound_id,
                    "smiles": smiles,
                    "corporate_id": corporate_id,
                }
            )
    return matches


def similarity_search(query_smiles: str, threshold: float = 0.7, db_path: Optional[Path] = None) -> List[dict]:
    """Similarity search using Tanimoto and RDKit fingerprints."""
    if not RDKIT_AVAILABLE:
        return []
    if not query_smiles:
        return []

    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        logger.warning("[CHEMISTRY][SEARCH] Invalid SMILES: %s", query_smiles)
        return []
    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)

    results: List[dict] = []
    for compound_id, smiles, corporate_id in _load_compounds(db_path):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        sim = DataStructs.TanimotoSimilarity(query_fp, fp)
        if sim >= threshold:
            results.append(
                {
                    "compound_id": compound_id,
                    "smiles": smiles,
                    "corporate_id": corporate_id,
                    "similarity": sim,
                }
            )

    # Sort by similarity descending
    results.sort(key=lambda x: x.get("similarity", 0), reverse=True)
    return results
