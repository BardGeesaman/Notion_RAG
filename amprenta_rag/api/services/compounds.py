from __future__ import annotations

from typing import Dict, List, Optional

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Compound


def list_compounds() -> List[Dict]:
    """List compounds ordered by latest update."""
    with db_session() as db:
        rows = (
            db.query(Compound)
            .order_by(Compound.updated_at.desc() if hasattr(Compound, "updated_at") else Compound.created_at.desc())
            .all()
        )
        return [_compound_to_dict(c) for c in rows]


def get_compound_by_id(compound_id: str) -> Optional[Dict]:
    """Fetch a compound by compound_id."""
    with db_session() as db:
        obj = db.query(Compound).filter(Compound.compound_id == compound_id).first()
        return _compound_to_dict(obj) if obj else None


def get_compound_programs(compound_id: str) -> List[Dict]:
    """Return programs linked to a compound."""
    with db_session() as db:
        obj = db.query(Compound).filter(Compound.compound_id == compound_id).first()
        if not obj or not getattr(obj, "programs", None):
            return []
        return [
            {
                "id": str(getattr(p, "id", "")) if getattr(p, "id", None) else None,
                "name": getattr(p, "name", None),
                "description": getattr(p, "description", None),
            }
            for p in obj.programs
        ]


def _compound_to_dict(c: Compound) -> Dict:
    if not c:
        return {}
    created_at_val = getattr(c, "created_at", None)
    updated_at_val = getattr(c, "updated_at", None)
    return {
        "id": str(getattr(c, "id", "")) if getattr(c, "id", None) else None,
        "compound_id": getattr(c, "compound_id", None),
        "smiles": getattr(c, "smiles", None),
        "inchi_key": getattr(c, "inchi_key", None),
        "canonical_smiles": getattr(c, "canonical_smiles", None),
        "molecular_formula": getattr(c, "molecular_formula", None),
        "molecular_weight": getattr(c, "molecular_weight", None),
        "logp": getattr(c, "logp", None),
        "hbd_count": getattr(c, "hbd_count", None),
        "hba_count": getattr(c, "hba_count", None),
        "rotatable_bonds": getattr(c, "rotatable_bonds", None),
        "created_at": created_at_val.isoformat() if created_at_val else None,
        "updated_at": updated_at_val.isoformat() if updated_at_val else None,
    }

