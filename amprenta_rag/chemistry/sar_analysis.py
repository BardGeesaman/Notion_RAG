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


def detect_activity_cliffs(threshold: float = 0.3) -> List[Dict[str, Any]]:
    """Detect activity cliffs (similar structures, different activities)."""
    # Placeholder - would need activity data
    return []
