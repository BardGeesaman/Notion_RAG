"""Seed synthetic metabolomics data with sparse features and adducts."""

from __future__ import annotations

import argparse
import random
from typing import Dict, List, Tuple

import numpy as np

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Feature


SIZE_PRESETS: Dict[str, int] = {
    "small": 30,
    "medium": 150,
    "large": 500,
}

ADDUCTS = ["[M+H]+", "[M+Na]+", "[M+K]+", "[M-H]-", "[M+Cl]-"]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed metabolomics demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo metabolites first.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    with db_session() as db:
        deleted = db.query(Feature).filter(Feature.name.like("MET_DEMO_%")).delete(synchronize_session=False)
        db.commit()
    return deleted


def _sparse_matrix(rows: int, cols: int, zero_fraction: float = 0.5) -> np.ndarray:
    mat = np.random.lognormal(mean=1.5, sigma=0.8, size=(rows, cols))
    mask = np.random.rand(rows, cols) < zero_fraction
    mat[mask] = 0.0
    return mat


def _mock_kegg_chebi(idx: int) -> Dict[str, str]:
    return {
        "kegg": f"C{idx:05d}",
        "chebi": f"CHEBI:{100000 + idx}",
    }


def _make_names(n: int) -> List[str]:
    return [f"MET_DEMO_{i+1:04d}" for i in range(n)]


def _seed_features(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int]:
    metabolites = SIZE_PRESETS[size]
    samples = 12 if size == "small" else 24 if size == "medium" else 48

    random.seed(seed_value)
    np.random.seed(seed_value)

    names = _make_names(metabolites)
    matrix = _sparse_matrix(metabolites, samples, zero_fraction=random.uniform(0.4, 0.6))
    mean_intensity = matrix.mean(axis=1)
    zeros = (matrix == 0).sum(axis=1)

    with db_session() as db:
        datasets = db.query(Dataset).filter(Dataset.name.like("DEMO_DS_%")).all()
        if not datasets:
            print("No DEMO datasets found. Run seed_core_scaffolding.py first.")
            return 0, 0

        created = linked = 0
        for i, name in enumerate(names):
            feature = db.query(Feature).filter(Feature.name == name).first()
            adduct = random.choice(ADDUCTS)
            mappings = _mock_kegg_chebi(i + 1)
            if not feature:
                feature = Feature(
                    name=name,
                    feature_type="metabolite",
                    normalized_name=name.lower(),
                    external_ids={
                        "adduct": adduct,
                        "mean_intensity": float(mean_intensity[i]),
                        "zero_count": int(zeros[i]),
                        "kegg_id": mappings["kegg"],
                        "chebi_id": mappings["chebi"],
                        "sample_count": samples,
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
        print(f"Reset: deleted demo metabolites = {deleted}")

    created, linked = _seed_features(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Features created: {created}")
    print(f"Feature-dataset links added: {linked}")


if __name__ == "__main__":
    main()

