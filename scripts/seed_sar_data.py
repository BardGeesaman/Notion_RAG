#!/usr/bin/env python3
"""
Seed deterministic SAR demo data (idempotent).

Creates:
- A small CDK2 SAR series (Compounds)
- BiochemicalResult rows with IC50 (nM)
- A guaranteed activity cliff pair (high similarity + large fold-change)

Usage:
  python scripts/seed_sar_data.py
  python scripts/seed_sar_data.py --reset
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

# Ensure repo root is on sys.path when running as `python scripts/seed_sar_data.py`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import BiochemicalResult, Compound


TARGET = "CDK2"
ASSAY_NAME = "CDK2 IC50 Assay"
UNITS = "nM"

# Stable prefixes so we can idempotently upsert/delete.
COMPOUND_PREFIX = "SAR-CDK2-"
RESULT_SUFFIX = "-RESULT"


def _now_utc() -> datetime:
    return datetime.now(timezone.utc)


def _inchi_key_best_effort(smiles: str) -> Optional[str]:
    """Return InChIKey if RDKit is available; otherwise None (best-effort parse)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
    except Exception:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        if mol is None:
            return None
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def _get_or_create_compound(
    *,
    db,
    compound_id: str,
    smiles: str,
    inchi_key: Optional[str],
) -> Compound:
    # Prefer exact compound_id match
    obj = db.query(Compound).filter(Compound.compound_id == compound_id).first()
    if obj:
        # Keep it simple: update SMILES if missing (don’t overwrite existing)
        if not getattr(obj, "smiles", None):
            obj.smiles = smiles
        if inchi_key and not getattr(obj, "inchi_key", None):
            obj.inchi_key = inchi_key
        return obj

    # Avoid inchi_key unique collisions by reusing an existing compound if present
    if inchi_key:
        existing = db.query(Compound).filter(Compound.inchi_key == inchi_key).first()
        if existing:
            return existing

    obj = Compound(
        compound_id=compound_id,
        smiles=smiles,
        inchi_key=inchi_key,
        canonical_smiles=None,
        created_at=_now_utc(),
        updated_at=_now_utc(),
    )
    db.add(obj)
    db.flush()
    return obj


def _upsert_biochemical_result(
    *,
    db,
    result_id: str,
    compound_uuid,
    target: str,
    assay_name: str,
    ic50_nm: float,
    units: str,
) -> BiochemicalResult:
    obj = db.query(BiochemicalResult).filter(BiochemicalResult.result_id == result_id).first()
    if obj:
        obj.compound_id = compound_uuid
        obj.target = target
        obj.assay_name = assay_name
        obj.ic50 = float(ic50_nm)
        obj.units = units
        obj.updated_at = _now_utc()
        return obj

    obj = BiochemicalResult(
        result_id=result_id,
        compound_id=compound_uuid,
        assay_name=assay_name,
        target=target,
        ic50=float(ic50_nm),
        units=units,
        run_date=_now_utc(),
        created_at=_now_utc(),
        updated_at=_now_utc(),
    )
    db.add(obj)
    return obj


def _seed_series(size: str) -> List[Tuple[str, str, float]]:
    """
    Return (compound_id, smiles, ic50_nm) rows.

    Notes:
    - SMILES are intentionally simple/aromatic for SAR demos.
    - One pair is a guaranteed activity cliff (similar enough, huge fold).
    """
    rows: List[Tuple[str, str, float]] = []

    # Deterministic cliff pair (halogen swap) with large delta
    rows.append((f"{COMPOUND_PREFIX}001", "Fc1ccc2[nH]c3ccccc3n2c1", 50.0))
    rows.append((f"{COMPOUND_PREFIX}002", "Clc1ccc2[nH]c3ccccc3n2c1", 1000.0))

    # Core heteroaromatic motif variants (baseline set)
    base_extra: List[Tuple[str, str, float]] = [
        ("COc1ccc2[nH]c3ccccc3n2c1", 120.0),
        ("Cc1ccc2[nH]c3ccccc3n2c1", 90.0),
        ("Nc1ccc2[nH]c3ccccc3n2c1", 300.0),
        ("Oc1ccc2[nH]c3ccccc3n2c1", 180.0),
        ("c1ccc2[nH]c3ccccc3n2c1", 250.0),
        ("Brc1ccc2[nH]c3ccccc3n2c1", 700.0),
    ]

    # Diverse scaffolds for larger sizes
    diverse: List[Tuple[str, str, float]] = [
        ("COc1nccc2ccccc12", 80.0),
        ("COc1nc(C)cc2ccccc12", 150.0),
        ("COc1nc(Cl)cc2ccccc12", 600.0),
        ("Cc1nccc2ccccc12", 110.0),
        ("Clc1nc(C)cc2ccccc12", 900.0),
        ("NC(=O)c1ccc2nc(C)ccc2c1", 400.0),
        ("COc1ccc2nccc(N)c2c1", 220.0),
        ("COc1ccc2nc(C)cc(Cl)c2c1", 750.0),
    ]

    if size == "small":
        extra = base_extra[:4]
    elif size == "medium":
        extra = base_extra + diverse[:4]
    else:
        extra = base_extra + diverse

    for idx, (smi, ic50) in enumerate(extra, start=3):
        rows.append((f"{COMPOUND_PREFIX}{idx:03d}", smi, float(ic50)))

    return rows


def seed(reset: bool = False, size: str = "small", dry_run: bool = False) -> Dict[str, int]:
    with db_session() as db:
        if reset:
            # Delete biochemical results first (FK points to compounds)
            deleted_results = (
                db.query(BiochemicalResult)
                .filter(BiochemicalResult.result_id.like(f"{COMPOUND_PREFIX}%{RESULT_SUFFIX}"))
                .delete(synchronize_session=False)
            )
            # Delete compounds (safe: unique compound_id prefix for seeded data)
            deleted_compounds = (
                db.query(Compound)
                .filter(Compound.compound_id.like(f"{COMPOUND_PREFIX}%"))
                .delete(synchronize_session=False)
            )
            db.commit()
        else:
            deleted_results = 0
            deleted_compounds = 0

        created_or_updated_compounds = 0
        created_or_updated_results = 0

        for compound_id, smiles, ic50 in _seed_series(size):
            inchi_key = _inchi_key_best_effort(smiles)
            c = _get_or_create_compound(db=db, compound_id=compound_id, smiles=smiles, inchi_key=inchi_key)
            created_or_updated_compounds += 1

            result_id = f"{compound_id}{RESULT_SUFFIX}"
            _upsert_biochemical_result(
                db=db,
                result_id=result_id,
                compound_uuid=c.id,
                target=TARGET,
                assay_name=ASSAY_NAME,
                ic50_nm=ic50,
                units=UNITS,
            )
            created_or_updated_results += 1

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return {
        "deleted_results": deleted_results,
        "deleted_compounds": deleted_compounds,
        "upserted_compounds": created_or_updated_compounds,
        "upserted_results": created_or_updated_results,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reset", action="store_true", help="Delete existing seeded SAR-CDK2-* rows first")
    ap.add_argument("--size", choices=["small", "medium", "large"], default="small", help="Size preset.")
    ap.add_argument("--seed", type=int, default=101, help="Random seed.")
    ap.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    args = ap.parse_args()

    # Re-seed PRNGs inside _seed_series via manual override if needed
    out = seed(reset=args.reset, size=args.size, dry_run=args.dry_run)
    print(f"SAR seed complete (size={args.size}, seed={args.seed}, dry_run={args.dry_run}):", out)


if __name__ == "__main__":
    main()


