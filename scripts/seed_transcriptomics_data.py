"""Seed synthetic transcriptomics data (genes, counts, differential expression)."""

from __future__ import annotations

import argparse
import math
import random
from typing import Dict, List, Tuple

import numpy as np

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset, Feature


SIZE_PRESETS: Dict[str, Tuple[int, int]] = {
    "small": (100, 6),    # genes, samples
    "medium": (500, 20),
    "large": (2000, 50),
}

COMMON_GENES = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "NRAS", "PIK3CA",
    "PTEN", "ALK", "BRAF", "ERBB2", "MET", "MYC", "CDKN2A", "RB1",
]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed transcriptomics demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo genes first.")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    with db_session() as db:
        deleted = db.query(Feature).filter(Feature.name.like("GENE_DEMO_%")).delete(synchronize_session=False)
        db.commit()
    return deleted


def _nbinom_counts(genes: int, samples: int) -> np.ndarray:
    # Negative binomial with moderate dispersion
    mean = 50
    dispersion = 0.4
    p = mean / (mean + dispersion)
    r = dispersion
    return np.random.negative_binomial(r, p, size=(genes, samples))


def _make_gene_names(n: int) -> List[str]:
    names: List[str] = []
    idx = 0
    # start with common genes
    while idx < n and idx < len(COMMON_GENES):
        names.append(COMMON_GENES[idx])
        idx += 1
    while idx < n:
        names.append(f"GENE_DEMO_{idx+1:05d}")
        idx += 1
    return names


def _fdr_from_pvalues(pvals: List[float]) -> List[float]:
    n = len(pvals)
    sorted_pairs = sorted(enumerate(pvals), key=lambda x: x[1])
    fdr = [0.0] * n
    for rank, (idx, p) in enumerate(sorted_pairs, start=1):
        fdr[idx] = min(1.0, p * n / rank)
    return fdr


def _seed_features(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int]:
    genes, samples = SIZE_PRESETS[size]
    random.seed(seed_value)
    np.random.seed(seed_value)

    gene_names = _make_gene_names(genes)
    counts = _nbinom_counts(genes, samples)
    log2fc = np.random.normal(loc=0.0, scale=1.2, size=genes)
    pvals = np.clip(np.random.uniform(low=0.0001, high=0.1, size=genes), 0, 1)
    fdrs = _fdr_from_pvalues(pvals.tolist())
    means = counts.mean(axis=1)

    # target datasets (created by scaffolding)
    with db_session() as db:
        datasets = db.query(Dataset).filter(Dataset.name.like("DEMO_DS_%")).all()
        if not datasets:
            print("No DEMO datasets found. Run seed_core_scaffolding.py first.")
            return 0, 0

        created = linked = 0
        for i, gene in enumerate(gene_names):
            feature = db.query(Feature).filter(Feature.name == gene).first()
            if not feature:
                feature = Feature(
                    name=gene,
                    feature_type="gene",
                    normalized_name=gene.lower(),
                    external_ids={
                        "log2fc": float(log2fc[i]),
                        "pvalue": float(pvals[i]),
                        "fdr": float(fdrs[i]),
                        "mean_count": float(means[i]),
                        "sample_count": samples,
                    },
                )
                db.add(feature)
                created += 1
            # link to all demo datasets
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
        print(f"Reset: deleted demo genes = {deleted}")

    created, linked = _seed_features(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"Features created: {created}")
    print(f"Feature-dataset links added: {linked}")


if __name__ == "__main__":
    main()

