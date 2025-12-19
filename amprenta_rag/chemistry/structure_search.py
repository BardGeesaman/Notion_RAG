"""Structure search utilities using PostgreSQL."""
from typing import List, Dict, Any

from amprenta_rag.logging_utils import get_logger

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound

logger = get_logger(__name__)

# Check RDKit availability
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def substructure_search(smarts: str, limit: int = 100) -> List[Dict[str, Any]]:
    """Search compounds by SMARTS substructure pattern using PostgreSQL."""
    if not RDKIT_AVAILABLE:
        return []

    try:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return []
    except Exception as e:
        logger.warning("[STRUCTURE_SEARCH] Failed to parse SMARTS '%s': %r", smarts, e)
        return []

    db_gen = get_db()
    db = next(db_gen)
    try:
        compounds = db.query(Compound).limit(1000).all()
        matches = []

        for c in compounds:
            if not c.smiles:
                continue
            mol = Chem.MolFromSmiles(c.smiles)
            if mol and mol.HasSubstructMatch(pattern):
                matches.append({
                    "compound_id": c.compound_id,
                    "smiles": c.smiles,
                    "molecular_weight": c.molecular_weight,
                })
                if len(matches) >= limit:
                    break

        return matches
    finally:
        db_gen.close()


def similarity_search(smiles: str, threshold: float = 0.7, limit: int = 100) -> List[Dict[str, Any]]:
    """Search compounds by Tanimoto similarity using PostgreSQL."""
    if not RDKIT_AVAILABLE:
        return []

    try:
        query_mol = Chem.MolFromSmiles(smiles)
        if query_mol is None:
            return []
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
    except Exception as e:
        logger.warning("[STRUCTURE_SEARCH] Failed to compute fingerprint for '%s': %r", smiles, e)
        return []

    db_gen = get_db()
    db = next(db_gen)
    try:
        compounds = db.query(Compound).limit(1000).all()
        results = []

        for c in compounds:
            if not c.smiles:
                continue
            mol = Chem.MolFromSmiles(c.smiles)
            if mol is None:
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(query_fp, fp)

            if similarity >= threshold:
                results.append({
                    "compound_id": c.compound_id,
                    "smiles": c.smiles,
                    "similarity": similarity,
                    "molecular_weight": c.molecular_weight,
                })

        results.sort(key=lambda x: x["similarity"], reverse=True)
        return results[:limit]
    finally:
        db_gen.close()
