"""SAR analysis utilities using PostgreSQL."""
from typing import Optional, List, Dict, Any

import pandas as pd

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound


def get_compound_properties(limit: int = 1000) -> pd.DataFrame:
    """Get compound properties from PostgreSQL."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        compounds = db.query(Compound).limit(limit).all()
        if not compounds:
            return pd.DataFrame()

        data = [{
            "compound_id": c.compound_id,
            "smiles": c.smiles,
            "molecular_weight": c.molecular_weight,
            "logp": c.logp,
            "hbd_count": c.hbd_count,
            "hba_count": c.hba_count,
            "rotatable_bonds": c.rotatable_bonds,
        } for c in compounds]
        return pd.DataFrame(data)
    finally:
        db_gen.close()


def get_activity_data(compound_ids: Optional[List[str]] = None) -> pd.DataFrame:
    """Get activity data for compounds from PostgreSQL."""
    db_gen = get_db()
    db = next(db_gen)
    try:
        query = db.query(Compound)
        if compound_ids:
            query = query.filter(Compound.compound_id.in_(compound_ids))
        compounds = query.all()

        if not compounds:
            return pd.DataFrame()

        data = [{
            "compound_id": c.compound_id,
            "smiles": c.smiles,
            "molecular_weight": c.molecular_weight,
        } for c in compounds]
        return pd.DataFrame(data)
    finally:
        db_gen.close()


def calculate_lipinski(smiles: str) -> Dict[str, Any]:
    """Calculate Lipinski Rule of 5 properties."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"valid": False}

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ])

        return {
            "valid": True,
            "molecular_weight": mw,
            "logp": logp,
            "hbd": hbd,
            "hba": hba,
            "violations": violations,
            "passes_ro5": violations == 0,
        }
    except ImportError:
        return {"valid": False, "error": "RDKit not available"}


def detect_activity_cliffs(
    assay_id: Optional[str] = None,
    similarity_threshold: float = 0.7,
    activity_fold_change: float = 10.0,
    limit: int = 50
) -> List[Dict[str, Any]]:
    """
    Detect activity cliffs - pairs of similar compounds with large activity differences.

    Args:
        assay_id: Optional assay UUID to filter results
        similarity_threshold: Minimum Tanimoto similarity (default 0.7)
        activity_fold_change: Minimum fold change in activity (default 10x)
        limit: Maximum number of cliffs to return

    Returns:
        List of activity cliff pairs with compound info and activity data
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except ImportError:
        return []

    from amprenta_rag.database.base import get_db
    from amprenta_rag.database.models import Compound, ActivityResult
    from uuid import UUID

    db_gen = get_db()
    db = next(db_gen)
    try:
        # Get compounds with activity data
        query = db.query(ActivityResult).join(Compound)
        if assay_id:
            query = query.filter(ActivityResult.assay_id == UUID(assay_id))

        results = query.all()
        if len(results) < 2:
            return []

        # Build compound data with fingerprints
        compound_data = {}
        for r in results:
            cid = str(r.compound_id)
            if cid not in compound_data:
                mol = Chem.MolFromSmiles(r.compound.smiles) if r.compound.smiles else None
                if mol:
                    compound_data[cid] = {
                        "compound_id": r.compound.compound_id,
                        "smiles": r.compound.smiles,
                        "fp": AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
                        "activities": {}
                    }
            if cid in compound_data:
                assay_key = str(r.assay_id)
                compound_data[cid]["activities"][assay_key] = r.value

        # Find activity cliffs
        cliffs = []
        compound_ids = list(compound_data.keys())

        for i, cid1 in enumerate(compound_ids):
            for cid2 in compound_ids[i+1:]:
                c1 = compound_data[cid1]
                c2 = compound_data[cid2]

                # Calculate similarity
                similarity = DataStructs.TanimotoSimilarity(c1["fp"], c2["fp"])

                if similarity >= similarity_threshold:
                    # Check for activity difference in shared assays
                    shared_assays = set(c1["activities"].keys()) & set(c2["activities"].keys())

                    for assay in shared_assays:
                        v1 = c1["activities"][assay]
                        v2 = c2["activities"][assay]

                        if v1 > 0 and v2 > 0:
                            fold_change = max(v1/v2, v2/v1)

                            if fold_change >= activity_fold_change:
                                cliffs.append({
                                    "compound_1": c1["compound_id"],
                                    "smiles_1": c1["smiles"],
                                    "activity_1": v1,
                                    "compound_2": c2["compound_id"],
                                    "smiles_2": c2["smiles"],
                                    "activity_2": v2,
                                    "similarity": similarity,
                                    "fold_change": fold_change,
                                    "assay_id": assay,
                                })

                                if len(cliffs) >= limit:
                                    return cliffs

        return cliffs
    finally:
        db_gen.close()
