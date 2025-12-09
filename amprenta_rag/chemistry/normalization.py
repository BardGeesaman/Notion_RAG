"""
Compound normalization using RDKit (if available) or fallback methods.

Normalizes SMILES, generates InChIKeys, and computes molecular descriptors.
"""

from __future__ import annotations

import hashlib
from typing import Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning(
        "[CHEMISTRY][NORM] RDKit not available. Using fallback normalization."
    )


def normalize_smiles(smiles: str) -> tuple[str, Optional[str], Optional[str]]:
    """
    Normalize a SMILES string.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Tuple of (canonical_smiles, inchi_key, molecular_formula)
        Returns (original_smiles, None, None) if RDKit is not available
    """
    if not smiles or not smiles.strip():
        return smiles, None, None
    
    if RDKIT_AVAILABLE:
        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol is None:
                logger.warning(
                    "[CHEMISTRY][NORM] Invalid SMILES: %s",
                    smiles,
                )
                return smiles, None, None
            
            # Generate canonical SMILES
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            
            # Generate InChI key
            inchi_key = None
            try:
                from rdkit.Chem import inchi
                inchi_string = inchi.MolToInchi(mol)
                if inchi_string:
                    # Generate InChI key from InChI string
                    inchi_key = inchi.InchiToInchiKey(inchi_string)
            except Exception as e:
                logger.debug(
                    "[CHEMISTRY][NORM] Could not generate InChI key: %r",
                    e,
                )
            
            # Generate molecular formula
            molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
            
            return canonical_smiles, inchi_key, molecular_formula
            
        except Exception as e:
            logger.warning(
                "[CHEMISTRY][NORM] Error normalizing SMILES %s: %r",
                smiles,
                e,
            )
            return smiles, None, None
    else:
        # Fallback: just return cleaned SMILES
        return smiles.strip(), None, None


def compute_molecular_descriptors(smiles: str) -> dict[str, Optional[float | int]]:
    """
    Compute molecular descriptors for a compound.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of descriptor values
    """
    descriptors = {
        "molecular_weight": None,
        "logp": None,
        "hbd_count": None,
        "hba_count": None,
        "rotatable_bonds": None,
    }
    
    # Early return for empty SMILES
    if not smiles or not smiles.strip():
        return descriptors
    
    if not RDKIT_AVAILABLE:
        return descriptors
    
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return descriptors
        
        descriptors["molecular_weight"] = Descriptors.MolWt(mol)
        descriptors["logp"] = Descriptors.MolLogP(mol)
        descriptors["hbd_count"] = Descriptors.NumHDonors(mol)
        descriptors["hba_count"] = Descriptors.NumHAcceptors(mol)
        descriptors["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        
    except Exception as e:
        logger.debug(
            "[CHEMISTRY][NORM] Error computing descriptors for %s: %r",
            smiles,
            e,
        )
    
    return descriptors


def generate_compound_id(smiles: str) -> str:
    """
    Generate a unique compound ID from SMILES.
    
    Uses canonical SMILES if RDKit is available, otherwise uses hash of input.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Unique compound ID
    """
    if RDKIT_AVAILABLE:
        canonical_smiles, _, _ = normalize_smiles(smiles)
        # Use hash of canonical SMILES as ID
        return hashlib.sha256(canonical_smiles.encode()).hexdigest()[:16]
    else:
        # Fallback: use hash of input SMILES
        return hashlib.sha256(smiles.strip().encode()).hexdigest()[:16]

