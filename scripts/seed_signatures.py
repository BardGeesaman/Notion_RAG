"""Seed synthetic cross-omics signatures."""

from __future__ import annotations

import argparse
import random
from typing import Dict, List

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Signature, SignatureComponent, Feature


SIZE_PRESETS: Dict[str, int] = {
    "small": 5,
    "medium": 20,
    "large": 50,
}

MODALITY_MAP = {
    "gene": "transcriptomics",
    "protein": "proteomics",
    "metabolite": "metabolomics",
    "lipid": "lipidomics",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed cross-omics signatures.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo signatures first.")
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    with db_session() as db:
        deleted_components = db.query(SignatureComponent).filter(SignatureComponent.feature_name.like("DEMO_%")).delete(synchronize_session=False)
        deleted = db.query(Signature).filter(Signature.name.like("DEMO_SIG_%")).delete(synchronize_session=False)
        db.commit()
    return deleted + deleted_components


def _pick_features(limit: int, rng: random.Random) -> List[Feature]:
    with db_session() as db:
        feats = (
            db.query(Feature)
            .filter(Feature.name.like("GENE_DEMO_%") | Feature.name.like("PROT_DEMO_%") | Feature.name.like("MET_DEMO_%") | Feature.name.like("LIP_DEMO_%"))
            .limit(limit * 5)
            .all()
        )
    rng.shuffle(feats)
    return feats[:limit]


def _create_signature(idx: int, components: List[Feature], rng: random.Random) -> Signature:
    name = f"DEMO_SIG_{idx:03d}"
    modalities = list({MODALITY_MAP.get(c.feature_type, c.feature_type) for c in components})
    sig = Signature(
        name=name,
        description=f"Synthetic cross-omics signature {idx}",
        modalities=modalities,
        created_by_id=None,
    )
    # build components
    for feat in components:
        direction = rng.choice(["up", "down"])
        weight = round(rng.uniform(0.5, 2.0), 2)
        conf = round(rng.uniform(0.6, 0.99), 2)
        sig.components.append(
            SignatureComponent(
                feature_id=feat.id,
                feature_name=feat.name,
                feature_type=feat.feature_type,
                direction=direction,
                weight=weight,
                confidence=conf if hasattr(SignatureComponent, "confidence") else None,  # tolerate older schema
            )
        )
    return sig


def _seed_signatures(size: str, seed_value: int, dry_run: bool) -> int:
    target = SIZE_PRESETS[size]
    rng = random.Random(seed_value)

    # choose per-signature feature counts
    if size == "small":
        comp_min, comp_max = 5, 8
    elif size == "medium":
        comp_min, comp_max = 8, 15
    else:
        comp_min, comp_max = 10, 20

    created = 0
    with db_session() as db:
        for idx in range(1, target + 1):
            name = f"DEMO_SIG_{idx:03d}"
            existing = db.query(Signature).filter(Signature.name == name).first()
            if existing:
                continue

            comp_count = rng.randint(comp_min, comp_max)
            feats = _pick_features(comp_count, rng)
            if not feats:
                print("No seeded features found (genes/proteins/metabolites/lipids). Run feature seed scripts first.")
                break

            sig = _create_signature(idx, feats, rng)
            db.add(sig)
            created += 1

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return created


def main() -> None:
    args = _parse_args()
    if args.reset:
        deleted = _reset_demo()
        print(f"Reset: deleted demo signatures/components = {deleted}")

    created = _seed_signatures(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Signatures created: {created}")


if __name__ == "__main__":
    main()

