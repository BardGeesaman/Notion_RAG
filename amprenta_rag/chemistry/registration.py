"""Compound registration and duplicate checking utilities."""
from __future__ import annotations

from typing import Optional

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Compound
from amprenta_rag.chemistry.normalization import normalize_smiles, compute_molecular_descriptors
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _get_next_corporate_id(db) -> str:
    """Generate next corporate ID (AMP-XXXXX) by querying existing compounds."""
    result = db.query(Compound).filter(
        Compound.compound_id.like("AMP-%")
    ).order_by(Compound.compound_id.desc()).first()
    
    if result:
        try:
            current_num = int(result.compound_id.split("-")[1])
            next_num = current_num + 1
        except (IndexError, ValueError):
            next_num = 1
    else:
        next_num = 1
    
    return f"AMP-{next_num:05d}"


def check_duplicate(smiles: str) -> Optional[str]:
    """Check for existing compound by SMILES/InChIKey. Returns compound_id if found."""
    if not smiles or not smiles.strip():
        return None

    canonical, inchi_key, _ = normalize_smiles(smiles)
    
    db_gen = get_db()
    db = next(db_gen)
    try:
        existing = db.query(Compound).filter(
            (Compound.canonical_smiles == canonical) |
            (Compound.smiles == canonical) |
            (Compound.inchi_key == inchi_key)
        ).first()
        return existing.compound_id if existing else None
    finally:
        db_gen.close()


def register_compound(
    name: str,
    smiles: str,
    salt_form: Optional[str] = None,
    batch_number: Optional[str] = None,
    parent_id: Optional[str] = None,
    registered_by: Optional[str] = None,
) -> Optional[str]:
    """Register a compound. Returns corporate_id or existing duplicate id."""
    if not smiles or not smiles.strip():
        return None

    existing = check_duplicate(smiles)
    if existing:
        logger.info("[CHEMISTRY][REG] Duplicate: %s", existing)
        return existing

    canonical, inchi_key, formula = normalize_smiles(smiles)
    descriptors = compute_molecular_descriptors(canonical or smiles)

    db_gen = get_db()
    db = next(db_gen)
    try:
        corporate_id = _get_next_corporate_id(db)
        
        compound = Compound(
            compound_id=corporate_id,
            smiles=canonical or smiles,
            inchi_key=inchi_key,
            canonical_smiles=canonical,
            molecular_formula=formula,
            molecular_weight=descriptors.get("molecular_weight"),
            logp=descriptors.get("logp"),
            hbd_count=descriptors.get("hbd_count"),
            hba_count=descriptors.get("hba_count"),
            rotatable_bonds=descriptors.get("rotatable_bonds"),
            aromatic_rings=descriptors.get("aromatic_rings"),
        )
        
        db.add(compound)
        db.commit()
        db.refresh(compound)
        logger.info("[CHEMISTRY][REG] Registered %s", corporate_id)
        
        # Fire workflow trigger
        from amprenta_rag.automation.engine import fire_trigger
        fire_trigger("compound_registered", {
            "compound_id": str(compound.id),
            "smiles": compound.smiles
        }, db)
        
        return corporate_id
    finally:
        db_gen.close()
