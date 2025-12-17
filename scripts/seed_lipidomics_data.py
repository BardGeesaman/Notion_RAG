"""Seed synthetic lipidomics data with chain-length/saturation annotations."""

from __future__ import annotations

import argparse
import random
from typing import Dict, List, Tuple

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Feature


SIZE_PRESETS: Dict[str, int] = {
    "small": 25,
    "medium": 100,
    "large": 400,
}

LIPID_CLASSES = ["Cer", "SM", "PC", "PE", "PI", "PS", "TG", "DG"]
CHAIN_BACKBONES = ["d18:1", "d18:0", "d20:1"]
FATTY_ACYLS = ["16:0", "18:0", "18:1", "20:4", "22:6", "24:0", "24:1"]
PATHWAYS = ["Sphingolipid metabolism", "Glycerophospholipid metabolism", "Fatty acid metabolism", "Triglyceride metabolism"]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed lipidomics demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo lipids first.")
    parser.add_argument("--seed", type=int, default=99)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    with db_session() as db:
        deleted = db.query(Feature).filter(Feature.name.like("LIP_DEMO_%")).delete(synchronize_session=False)
        db.commit()
    return deleted


def _make_lipid_names(n: int, rng: random.Random) -> List[Dict[str, str]]:
    lipids: List[Dict[str, str]] = []
    for i in range(n):
        cls = rng.choice(LIPID_CLASSES)
        backbone = rng.choice(CHAIN_BACKBONES)
        fa = rng.choice(FATTY_ACYLS)
        name = f"{cls}({backbone}/{fa})"
        lipids.append(
            {
                "name": f"LIP_DEMO_{i+1:04d}",
                "display_name": name,
                "class": cls,
                "backbone": backbone,
                "fa": fa,
                "pathway": rng.choice(PATHWAYS),
            }
        )
    return lipids


def _seed_features(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int]:
    total = SIZE_PRESETS[size]
    rng = random.Random(seed_value)

    lipid_defs = _make_lipid_names(total, rng)

    with db_session() as db:
        datasets = db.query(Dataset).filter(Dataset.name.like("DEMO_DS_%")).all()
        if not datasets:
            print("No DEMO datasets found. Run seed_core_scaffolding.py first.")
            return 0, 0

        created = linked = 0
        for lip in lipid_defs:
            feature = db.query(Feature).filter(Feature.name == lip["name"]).first()
            if not feature:
                feature = Feature(
                    name=lip["name"],
                    feature_type="lipid",
                    normalized_name=lip["display_name"].lower(),
                    external_ids={
                        "display_name": lip["display_name"],
                        "class": lip["class"],
                        "backbone": lip["backbone"],
                        "fatty_acyl": lip["fa"],
                        "pathway": lip["pathway"],
                    },
                )
                db.add(feature)
                created += 1
            for ds in datasets:
                if feature not in ds.features:
                    ds.features.append(feature)
                    linked += 1

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return created, linked


def main() -> None:
    args = _parse_args()
    if args.reset:
        deleted = _reset_demo()
        print(f"Reset: deleted demo lipids = {deleted}")

    created, linked = _seed_features(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Features created: {created}")
    print(f"Feature-dataset links added: {linked}")


if __name__ == "__main__":
    main()

