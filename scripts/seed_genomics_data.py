"""Seed synthetic genomics (variant) data (VariantSet, Variant, VariantAnnotation)."""

from __future__ import annotations

import argparse
import random
from typing import Dict, List, Tuple

from amprenta_rag.database.models import Variant, VariantAnnotation, VariantSet
from amprenta_rag.database.session import db_session


SIZE_PRESETS: Dict[str, Tuple[int, int, int]] = {
    # sets, variants, annotations
    "small": (1, 100, 20),
    "medium": (5, 1_000, 200),
    "large": (20, 10_000, 2_000),
}

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
BASES = ["A", "C", "G", "T"]
COMMON_GENES = [
    "TP53",
    "BRCA1",
    "BRCA2",
    "EGFR",
    "KRAS",
    "NRAS",
    "PIK3CA",
    "PTEN",
    "ALK",
    "BRAF",
    "ERBB2",
    "MET",
    "MYC",
    "CDKN2A",
    "RB1",
]

CONSEQUENCES = [
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_region_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
]
IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
CLIN_SIG = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed genomics (variant) demo data.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo variant sets first.")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _reset_demo() -> int:
    """Delete demo variant sets (cascades to variants + annotations)."""
    with db_session() as db:
        sets = db.query(VariantSet).filter(VariantSet.name.like("DEMO_VARSET_%")).all()
        deleted = len(sets)
        for vs in sets:
            db.delete(vs)
        db.commit()
    return deleted


def _make_gene_symbol(idx: int) -> str:
    if idx < len(COMMON_GENES):
        return COMMON_GENES[idx]
    return f"GENE_DEMO_{idx+1:05d}"


def _pick_alleles(rng: random.Random) -> Tuple[str, str]:
    ref = rng.choice(BASES)
    alt = rng.choice([b for b in BASES if b != ref])
    return ref, alt


def _hgvs_genomic(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"chr{chrom}:g.{pos}{ref}>{alt}"


def _seed(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int, int]:
    rng = random.Random(seed_value)
    n_sets, n_vars, n_ann = SIZE_PRESETS[size]

    created_sets = created_vars = created_anns = 0

    # Distribute variants approximately evenly across sets.
    base_per_set = n_vars // n_sets
    remainder = n_vars % n_sets

    # Distribute annotations across the first N variants we generate.
    ann_remaining = n_ann

    with db_session() as db:
        for set_idx in range(n_sets):
            this_n = base_per_set + (1 if set_idx < remainder else 0)
            name = f"DEMO_VARSET_{set_idx+1:03d}"
            existing = db.query(VariantSet).filter(VariantSet.name == name).first()
            if existing:
                # Idempotency: do not duplicate VariantSets or their children on re-run without --reset.
                continue

            vs = VariantSet(
                name=name,
                description=f"Demo VariantSet ({size}) #{set_idx+1}",
                source_file=f"demo_{size}_{set_idx+1:03d}.vcf",
                source_type="vcf",
                status="completed",
            )
            db.add(vs)
            db.flush()  # assign vs.id
            created_sets += 1

            seen_loci = set()
            gene_symbols: List[str] = []

            for i in range(this_n):
                chrom = rng.choice(CHROMS)
                pos = rng.randint(10_000, 200_000_000)
                ref, alt = _pick_alleles(rng)
                locus = (chrom, pos, ref, alt)
                # Avoid collisions with the unique constraint per set.
                attempts = 0
                while locus in seen_loci and attempts < 10:
                    pos = rng.randint(10_000, 200_000_000)
                    ref, alt = _pick_alleles(rng)
                    locus = (chrom, pos, ref, alt)
                    attempts += 1
                seen_loci.add(locus)

                gene = _make_gene_symbol(rng.randint(0, 200))
                gene_symbols.append(gene)

                consequence = rng.choice(CONSEQUENCES)
                impact = rng.choice(IMPACTS)
                af = 10 ** rng.uniform(-6, -1)  # ~1e-6 to 1e-1

                var = Variant(
                    variant_set_id=vs.id,
                    chromosome=chrom,
                    position=pos,
                    ref_allele=ref,
                    alt_allele=alt,
                    rs_id=f"rs{rng.randint(1_000_000, 99_999_999)}" if rng.random() < 0.4 else None,
                    hgvs_genomic=_hgvs_genomic(chrom, pos, ref, alt),
                    hgvs_coding=None,
                    hgvs_protein=None,
                    gene_symbol=gene,
                    consequence=consequence,
                    impact=impact,
                    gnomad_af=af,
                )
                db.add(var)
                db.flush()  # assign var.id for annotations
                created_vars += 1

                if ann_remaining > 0 and rng.random() < 0.25:
                    ann = VariantAnnotation(
                        variant_id=var.id,
                        annotation_source=rng.choice(["clinvar", "vep"]),
                        clinvar_id=f"VCV{rng.randint(100_000, 999_999)}" if rng.random() < 0.5 else None,
                        clinical_significance=rng.choice(CLIN_SIG),
                        review_status=rng.choice(
                            [
                                "criteria provided, single submitter",
                                "reviewed by expert panel",
                                "no assertion criteria provided",
                            ]
                        ),
                        condition=rng.choice(
                            [
                                "ALS",
                                "Parkinson disease",
                                "Alzheimer disease",
                                "Frontotemporal dementia",
                                None,
                            ]
                        ),
                        cadd_phred=rng.uniform(1.0, 35.0),
                        revel_score=rng.uniform(0.0, 1.0),
                        sift_prediction=rng.choice(["tolerated", "deleterious", None]),
                        polyphen_prediction=rng.choice(["benign", "possibly_damaging", "probably_damaging", None]),
                    )
                    db.add(ann)
                    created_anns += 1
                    ann_remaining -= 1

            # Set summary counts
            vs.n_variants = this_n
            vs.n_genes = len(set(gene_symbols))

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return created_sets, created_vars, created_anns


def main() -> None:
    args = _parse_args()
    if args.reset:
        deleted = _reset_demo()
        print(f"Reset: deleted demo VariantSets = {deleted}")

    created_sets, created_vars, created_anns = _seed(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"VariantSets created: {created_sets}")
    print(f"Variants created: {created_vars}")
    print(f"VariantAnnotations created: {created_anns}")


if __name__ == "__main__":
    main()


