"""Pharmacophore feature extraction and search utilities."""
from __future__ import annotations

from typing import List, Optional

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    Chem = None  # type: ignore
    Descriptors = None  # type: ignore


def get_pharmacophore_features(smiles: str) -> dict:
    """
    Extract pharmacophore features from a SMILES string.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dict with pharmacophore features:
        - hbd: Number of hydrogen bond donors
        - hba: Number of hydrogen bond acceptors
        - aromatic_rings: Number of aromatic rings
        - hydrophobic: Hydrophobicity score (based on LogP/TPSA)
    """
    if not RDKIT_AVAILABLE:
        logger.error("[PHARMACOPHORE] RDKit not available")
        return {
            "hbd": 0,
            "hba": 0,
            "aromatic_rings": 0,
            "hydrophobic": 0.0,
        }
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("[PHARMACOPHORE] Invalid SMILES: %s", smiles)
            return {
                "hbd": 0,
                "hba": 0,
                "aromatic_rings": 0,
                "hydrophobic": 0.0,
            }
        
        # Calculate descriptors
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        
        # Calculate hydrophobicity score
        # Higher LogP and lower TPSA = more hydrophobic
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Normalize hydrophobicity: logP contributes positively, TPSA negatively
        # Scale TPSA to similar range as LogP (divide by 10)
        hydrophobic_score = logp - (tpsa / 10.0)
        
        features = {
            "hbd": int(hbd),
            "hba": int(hba),
            "aromatic_rings": int(aromatic_rings),
            "hydrophobic": float(hydrophobic_score),
        }
        
        logger.debug("[PHARMACOPHORE] Features for %s: %s", smiles[:30], features)
        return features
    
    except Exception as e:
        logger.error("[PHARMACOPHORE] Error extracting features: %r", e)
        return {
            "hbd": 0,
            "hba": 0,
            "aromatic_rings": 0,
            "hydrophobic": 0.0,
        }


def pharmacophore_search(
    min_hbd: Optional[int] = None,
    max_hbd: Optional[int] = None,
    min_hba: Optional[int] = None,
    max_hba: Optional[int] = None,
    min_aromatic: Optional[int] = None,
    max_aromatic: Optional[int] = None,
    min_hydrophobic: Optional[float] = None,
    max_hydrophobic: Optional[float] = None,
    db=None,
    limit: int = 100,
) -> List[Compound]:
    """
    Search compounds by pharmacophore features.
    
    Args:
        min_hbd: Minimum hydrogen bond donors
        max_hbd: Maximum hydrogen bond donors
        min_hba: Minimum hydrogen bond acceptors
        max_hba: Maximum hydrogen bond acceptors
        min_aromatic: Minimum aromatic rings
        max_aromatic: Maximum aromatic rings
        min_hydrophobic: Minimum hydrophobicity score
        max_hydrophobic: Maximum hydrophobicity score
        db: Database session (if None, creates new session)
        limit: Maximum number of results
        
    Returns:
        List of Compound objects matching the criteria
    """
    if not RDKIT_AVAILABLE:
        logger.error("[PHARMACOPHORE] RDKit not available")
        return []
    
    use_external_db = db is None
    if use_external_db:
        db_gen = get_db()
        db = next(db_gen)
    
    try:
        query = db.query(Compound)
        
        # Filter by HBD
        if min_hbd is not None:
            query = query.filter(Compound.hbd_count >= min_hbd)
        if max_hbd is not None:
            query = query.filter(Compound.hbd_count <= max_hbd)
        
        # Filter by HBA
        if min_hba is not None:
            query = query.filter(Compound.hba_count >= min_hba)
        if max_hba is not None:
            query = query.filter(Compound.hba_count <= max_hba)
        
        # Filter by aromatic rings
        # Note: We don't have aromatic_rings stored directly, so we'll need to calculate
        # For now, we'll skip this filter or approximate using other features
        # TODO: Add aromatic_rings column to Compound model or calculate on-the-fly
        
        # Filter by hydrophobicity (using LogP as proxy)
        if min_hydrophobic is not None:
            # Approximate: hydrophobic_score ≈ logp - (tpsa/10)
            # For simplicity, use logp as proxy
            query = query.filter(Compound.logp >= min_hydrophobic)
        if max_hydrophobic is not None:
            query = query.filter(Compound.logp <= max_hydrophobic)
        
        compounds = query.limit(limit).all()
        
        # Post-filter by aromatic rings if needed (requires RDKit calculation)
        if min_aromatic is not None or max_aromatic is not None:
            filtered = []
            for compound in compounds:
                if not compound.smiles:
                    continue
                features = get_pharmacophore_features(compound.smiles)
                aromatic_count = features["aromatic_rings"]
                
                if min_aromatic is not None and aromatic_count < min_aromatic:
                    continue
                if max_aromatic is not None and aromatic_count > max_aromatic:
                    continue
                
                filtered.append(compound)
            compounds = filtered
        
        logger.info("[PHARMACOPHORE] Found %d compounds matching criteria", len(compounds))
        return compounds
    
    except Exception as e:
        logger.error("[PHARMACOPHORE] Error searching compounds: %r", e)
        return []
    
    finally:
        if use_external_db:
            db_gen.close()
