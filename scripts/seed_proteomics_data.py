"""Seed synthetic proteomics data (protein intensities with missingness)."""

from __future__ import annotations

import argparse
import random
from typing import Dict, List, Tuple

import numpy as np

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Feature


SIZE_PRESETS: Dict[str, Tuple[int, int]] = {
    "small": (50, 6),    # proteins, samples
    "medium": (200, 20),
    "large": (1000, 50),
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed proteomics demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo proteins first.")
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    with db_session() as db:
        deleted = db.query(Feature).filter(Feature.name.like("PROT_DEMO_%")).delete(synchronize_session=False)
        db.commit()
    return deleted


def _lognormal_intensities(proteins: int, samples: int) -> np.ndarray:
    # log-normal with moderate variance
    mu = 2.0
    sigma = 0.7
    return np.random.lognormal(mean=mu, sigma=sigma, size=(proteins, samples))


def _apply_missingness(matrix: np.ndarray, missing_fraction: float) -> np.ndarray:
    mask = np.random.rand(*matrix.shape) < missing_fraction
    matrix_missing = matrix.copy()
    matrix_missing[mask] = np.nan
    return matrix_missing


def _make_protein_names(n: int) -> List[str]:
    return [f"PROT_DEMO_{i+1:05d}" for i in range(n)]


def _seed_features(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int]:
    proteins, samples = SIZE_PRESETS[size]
    random.seed(seed_value)
    np.random.seed(seed_value)

    names = _make_protein_names(proteins)
    intensities = _lognormal_intensities(proteins, samples)
    missing_fraction = random.uniform(0.1, 0.3)
    intensities = _apply_missingness(intensities, missing_fraction)
    means = np.nanmean(intensities, axis=1)
    missing_counts = np.isnan(intensities).sum(axis=1)

    with db_session() as db:
        datasets = db.query(Dataset).filter(Dataset.name.like("DEMO_DS_%")).all()
        if not datasets:
            print("No DEMO datasets found. Run seed_core_scaffolding.py first.")
            return 0, 0

        created = linked = 0
        for i, name in enumerate(names):
            feature = db.query(Feature).filter(Feature.name == name).first()
            if not feature:
                feature = Feature(
                    name=name,
                    feature_type="protein",
                    normalized_name=name.lower(),
                    external_ids={
                        "mean_intensity": float(means[i]),
                        "missing_values": int(missing_counts[i]),
                        "sample_count": samples,
                        "missing_fraction": float(missing_fraction),
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
        print(f"Reset: deleted demo proteins = {deleted}")

    created, linked = _seed_features(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Features created: {created}")
    print(f"Feature-dataset links added: {linked}")


if __name__ == "__main__":
    main()

