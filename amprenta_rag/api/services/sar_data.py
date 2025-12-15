from __future__ import annotations

from typing import Any, Dict, List

from sqlalchemy import desc, distinct, func

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import BiochemicalResult, Compound


def _mol_from_smiles_best_effort(smiles: str):
    """
    RDKit parser that tolerates "can't kekulize" demo SMILES.

    Some seeded aromatic systems can fail full sanitization; for fingerprints we can
    still parse with `sanitize=False` and run a partial sanitize that skips kekulization.
    """
    if not smiles:
        return None
    try:
        from rdkit import Chem
    except Exception:
        return None

    m = Chem.MolFromSmiles(smiles)
    if m is not None:
        return m

    # Fallback: tolerate kekulization issues
    try:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return None
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        return m
    except Exception:
        return None


def list_targets(limit: int = 200) -> List[Dict[str, Any]]:
    """
    List SAR targets derived from biochemical results.

    Returns objects shaped like:
      {"target": "CDK2", "compound_count": 42}
    """
    with db_session() as db:
        # Group by target; count distinct compounds per target
        rows = (
            db.query(
                BiochemicalResult.target,
                func.count(distinct(BiochemicalResult.compound_id)).label("compound_count"),
            )
            .filter(BiochemicalResult.target.isnot(None))
            .filter(BiochemicalResult.target != "")
            .group_by(BiochemicalResult.target)
            .order_by(desc("compound_count"))
            .limit(limit)
            .all()
        )

        out: List[Dict[str, Any]] = []
        for target, compound_count in rows:
            out.append({"target": target, "compound_count": int(compound_count or 0)})
        return out


def get_compounds_by_target(target: str, limit: int = 2000) -> List[Dict[str, Any]]:
    """
    Return compounds + activity values for a target.

    For now, uses `BiochemicalResult.ic50` as the activity value (seed data uses IC50 nM).
    """
    with db_session() as db:
        rows = (
            db.query(BiochemicalResult, Compound)
            .join(Compound, BiochemicalResult.compound_id == Compound.id)
            .filter(BiochemicalResult.target == target)
            .order_by(BiochemicalResult.ic50.asc().nullslast())
            .limit(limit)
            .all()
        )

        out: List[Dict[str, Any]] = []
        for br, c in rows:
            out.append(
                {
                    "compound_id": getattr(c, "compound_id", None),
                    "smiles": getattr(c, "smiles", None),
                    "ic50": getattr(br, "ic50", None),
                    "units": getattr(br, "units", None),
                    "assay_name": getattr(br, "assay_name", None),
                    "result_id": getattr(br, "result_id", None),
                }
            )
        # Filter out rows with missing essentials
        return [r for r in out if r.get("compound_id") and r.get("smiles")]


def get_activity_cliffs_for_target(
    target: str,
    similarity_threshold: float = 0.6,
    fold_change: float = 10.0,
    limit: int = 50,
) -> List[Dict[str, Any]]:
    """
    Detect activity cliffs within a target using Morgan fingerprints (RDKit).

    Uses IC50 values from `BiochemicalResult.ic50`. Returns a list of pairs:
      {"compound_1", "smiles_1", "activity_1", "compound_2", "smiles_2", "activity_2",
       "similarity", "fold_change", "assay_id": None}
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except Exception:
        return []

    rows = get_compounds_by_target(target, limit=5000)
    # Keep only rows with numeric, positive activity
    data = []
    for r in rows:
        v = r.get("ic50")
        try:
            fv = float(v) if v is not None else None
        except Exception:
            fv = None
        if fv is None or fv <= 0:
            continue
        data.append({**r, "ic50": fv})

    if len(data) < 2:
        return []

    fp: Dict[str, Any] = {}
    for r in data:
        smi = r.get("smiles")
        m = _mol_from_smiles_best_effort(smi)
        if not m:
            continue
        fp[r["compound_id"]] = (r, AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))

    ids = list(fp.keys())
    cliffs: List[Dict[str, Any]] = []

    for i, id1 in enumerate(ids):
        r1, fp1 = fp[id1]
        for id2 in ids[i + 1 :]:
            r2, fp2 = fp[id2]
            sim = float(DataStructs.TanimotoSimilarity(fp1, fp2))
            if sim < similarity_threshold:
                continue
            v1 = float(r1["ic50"])
            v2 = float(r2["ic50"])
            fc = max(v1 / v2, v2 / v1)
            if fc < fold_change:
                continue
            cliffs.append(
                {
                    "compound_1": r1["compound_id"],
                    "smiles_1": r1["smiles"],
                    "activity_1": v1,
                    "compound_2": r2["compound_id"],
                    "smiles_2": r2["smiles"],
                    "activity_2": v2,
                    "similarity": round(sim, 3),
                    "fold_change": round(fc, 2),
                    "assay_id": None,
                }
            )
            if len(cliffs) >= limit:
                return cliffs

    # Sort best-first for stable UX
    cliffs.sort(key=lambda x: (x.get("fold_change", 0), x.get("similarity", 0)), reverse=True)
    return cliffs[:limit]


