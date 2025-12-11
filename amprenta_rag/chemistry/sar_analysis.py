"""Structure-Activity Relationship (SAR) analysis utilities."""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from amprenta_rag.chemistry.database import get_chemistry_db, get_chemistry_db_path
from amprenta_rag.chemistry.structure_search import RDKIT_AVAILABLE
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _get_db_path(db_path: Optional[Path] = None) -> Path:
    return db_path or get_chemistry_db_path()


def get_compound_properties(db_path: Optional[Path] = None) -> pd.DataFrame:
    """Load compound properties for SAR analysis."""
    path = _get_db_path(db_path)
    conn = get_chemistry_db(path)
    try:
        df = pd.read_sql(
            """
            SELECT compound_id, corporate_id, smiles, molecular_weight, logp, hbd_count, hba_count, rotatable_bonds
            FROM compounds
            """,
            conn,
        )
        return df
    finally:
        conn.close()


def get_activity_data(campaign_id: Optional[str] = None, db_path: Optional[Path] = None) -> pd.DataFrame:
    """Join compounds with HTS results to get activity data."""
    path = _get_db_path(db_path)
    conn = get_chemistry_db(path)
    try:
        base_sql = """
            SELECT c.compound_id, c.corporate_id, h.raw_value as activity_value, h.hit_flag
            FROM compounds c
            JOIN hts_results h ON c.compound_id = h.compound_id
        """
        params = ()
        if campaign_id:
            base_sql += " WHERE h.campaign_id = ?"
            params = (campaign_id,)
        df = pd.read_sql(base_sql, conn, params=params)
        return df
    finally:
        conn.close()


def calculate_lipinski(mw: Optional[float], logp: Optional[float], hbd: Optional[int], hba: Optional[int]) -> dict:
    """Evaluate Lipinski Rule of Five compliance."""
    violations = 0
    details = []
    if mw is not None and mw > 500:
        violations += 1
        details.append("MW > 500")
    if logp is not None and logp > 5:
        violations += 1
        details.append("LogP > 5")
    if hbd is not None and hbd > 5:
        violations += 1
        details.append("HBD > 5")
    if hba is not None and hba > 10:
        violations += 1
        details.append("HBA > 10")
    return {
        "compliant": violations == 0,
        "violations": violations,
        "details": details,
    }


def detect_activity_cliffs(activity_df: pd.DataFrame, similarity_threshold: float = 0.8) -> List[dict]:
    """Find activity cliffs: similar compounds with large activity differences."""
    if activity_df is None or activity_df.empty:
        return []
    if not RDKIT_AVAILABLE:
        logger.warning("[CHEMISTRY][SAR] RDKit not available; cannot detect activity cliffs.")
        return []

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
    except ImportError:
        return []

    # Precompute fingerprints
    fps = {}
    for _, row in activity_df.iterrows():
        smi = row.get("smiles") or ""
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fps[row["compound_id"]] = (AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048), row)

    cliffs: List[dict] = []
    ids = list(fps.keys())
    for i in range(len(ids)):
        fp_i, row_i = fps[ids[i]]
        for j in range(i + 1, len(ids)):
            fp_j, row_j = fps[ids[j]]
            sim = DataStructs.TanimotoSimilarity(fp_i, fp_j)
            if sim >= similarity_threshold:
                ai = row_i.get("activity_value")
                aj = row_j.get("activity_value")
                if ai is None or aj is None:
                    continue
                diff = abs(ai - aj)
                cliffs.append(
                    {
                        "compound_a": row_i.get("compound_id"),
                        "compound_b": row_j.get("compound_id"),
                        "similarity": sim,
                        "activity_diff": diff,
                    }
                )

    cliffs.sort(key=lambda x: x["activity_diff"], reverse=True)
    return cliffs
