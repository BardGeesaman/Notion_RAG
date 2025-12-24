"""Seed synthetic single-cell metadata (SingleCellDataset, CellCluster, CellTypeMarker)."""

from __future__ import annotations

import argparse
import random
from datetime import datetime, timezone
from typing import Dict, List, Tuple

from amprenta_rag.database.models import CellCluster, CellTypeMarker, Dataset, SingleCellDataset
from amprenta_rag.database.session import db_session


SIZE_PRESETS: Dict[str, Tuple[int, int, int, int]] = {
    # datasets, cells, clusters, markers
    "small": (2, 500, 5, 10),
    "medium": (10, 5_000, 20, 50),
    "large": (50, 50_000, 100, 200),
}

CELL_TYPES = [
    "T cell",
    "B cell",
    "NK cell",
    "Monocyte",
    "Macrophage",
    "Dendritic cell",
    "Endothelial",
    "Epithelial",
    "Fibroblast",
    "Astrocyte",
    "Neuron",
    "Microglia",
    "Oligodendrocyte",
]

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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Seed single-cell demo metadata.")
    parser.add_argument("--size", choices=SIZE_PRESETS.keys(), default="small")
    parser.add_argument("--reset", action="store_true", help="Delete existing demo single-cell datasets first.")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true", help="Simulate without committing changes.")
    return parser.parse_args()


def _distribute(total: int, parts: int) -> List[int]:
    if parts <= 0:
        return []
    base = total // parts
    rem = total % parts
    return [base + (1 if i < rem else 0) for i in range(parts)]


def _reset_demo() -> Tuple[int, int]:
    """Delete demo single-cell datasets and their backing Dataset rows."""
    with db_session() as db:
        # Find backing datasets
        demo_datasets = db.query(Dataset).filter(Dataset.name.like("DEMO_SC_DS_%")).all()
        ds_ids = [d.id for d in demo_datasets]

        # Delete single-cell datasets first (cascades clusters/markers/annotations)
        scds = []
        if ds_ids:
            scds = db.query(SingleCellDataset).filter(SingleCellDataset.dataset_id.in_(ds_ids)).all()
            for scd in scds:
                db.delete(scd)

        # Then delete datasets
        for d in demo_datasets:
            db.delete(d)

        db.commit()

    return len(demo_datasets), len(scds)


def _gene_symbol(rng: random.Random, idx: int) -> str:
    if idx < len(COMMON_GENES):
        return COMMON_GENES[idx]
    return f"SC_GENE_DEMO_{idx+1:05d}"


def _seed(size: str, seed_value: int, dry_run: bool) -> Tuple[int, int, int]:
    rng = random.Random(seed_value)
    n_datasets, total_cells, total_clusters, total_markers = SIZE_PRESETS[size]

    # Ensure at least one cluster per dataset
    base_clusters = [1] * n_datasets
    remaining_clusters = max(0, total_clusters - n_datasets)
    extra = _distribute(remaining_clusters, n_datasets)
    clusters_per_ds = [base_clusters[i] + extra[i] for i in range(n_datasets)]

    cells_per_ds = _distribute(total_cells, n_datasets)

    # Distribute markers roughly proportional to number of clusters in dataset
    markers_per_ds = [0] * n_datasets
    if total_markers > 0:
        weights = [max(1, c) for c in clusters_per_ds]
        wsum = sum(weights)
        assigned = 0
        for i in range(n_datasets):
            m = int(round(total_markers * (weights[i] / wsum)))
            markers_per_ds[i] = m
            assigned += m
        # Fix rounding drift
        while assigned > total_markers:
            j = rng.randrange(n_datasets)
            if markers_per_ds[j] > 0:
                markers_per_ds[j] -= 1
                assigned -= 1
        while assigned < total_markers:
            j = rng.randrange(n_datasets)
            markers_per_ds[j] += 1
            assigned += 1

    created_scd = created_clusters = created_markers = 0

    with db_session() as db:
        now = datetime.now(timezone.utc)
        gene_idx = 0
        for ds_idx in range(n_datasets):
            name = f"DEMO_SC_DS_{ds_idx+1:03d}"
            existing_ds = db.query(Dataset).filter(Dataset.name == name).first()
            if existing_ds:
                # Idempotency: do not duplicate Dataset/SingleCellDataset on re-run without --reset.
                continue

            ds = Dataset(
                name=name,
                omics_type="single_cell",
                description=f"Demo single-cell dataset ({size}) #{ds_idx+1}",
                ingestion_status="complete",
                external_ids={"demo": True, "seed": seed_value, "size": size},
            )
            db.add(ds)
            db.flush()

            n_cells = cells_per_ds[ds_idx]
            n_clusters = clusters_per_ds[ds_idx]
            n_genes = rng.randint(2000, 6000)
            hvg = min(3000, max(500, int(n_genes * 0.4)))

            scd = SingleCellDataset(
                dataset_id=ds.id,
                h5ad_path=f"data/single_cell/{ds.name}.h5ad",
                file_size_bytes=int(n_cells * n_genes * 4),
                n_cells=n_cells,
                n_genes=n_genes,
                normalization_method="log1p",
                hvg_count=hvg,
                clustering_method="leiden",
                clustering_resolution=0.5,
                has_pca=True,
                has_umap=True,
                processing_status="completed",
                processing_log="Seeded demo single-cell metadata.",
                ingested_at=now,
                processed_at=now,
            )
            db.add(scd)
            db.flush()
            created_scd += 1

            # Cluster cell counts should sum to dataset cells
            cells_per_cluster = _distribute(n_cells, n_clusters)
            for cid in range(n_clusters):
                cell_type = rng.choice(CELL_TYPES)
                cc = CellCluster(
                    single_cell_dataset_id=scd.id,
                    cluster_id=cid,
                    cell_type=cell_type,
                    n_cells=cells_per_cluster[cid],
                    marker_feature_ids=None,
                    description=f"Cluster {cid} ({cell_type})",
                )
                db.add(cc)
                created_clusters += 1

            # Markers
            mcount = markers_per_ds[ds_idx]
            for _ in range(mcount):
                cid = rng.randrange(n_clusters)
                gene = _gene_symbol(rng, gene_idx)
                gene_idx += 1

                pval = 10 ** rng.uniform(-6, -1)  # 1e-6 to 1e-1
                pval_adj = min(1.0, pval * 10.0)
                pct_in = rng.uniform(0.1, 1.0)
                pct_out = rng.uniform(0.0, min(0.95, pct_in))

                mk = CellTypeMarker(
                    single_cell_dataset_id=scd.id,
                    cluster_id=cid,
                    feature_id=None,
                    gene_symbol=gene,
                    log2_fold_change=rng.normalvariate(1.5, 0.75),
                    pval=pval,
                    pval_adj=pval_adj,
                    pct_in_cluster=pct_in,
                    pct_out_cluster=pct_out,
                )
                db.add(mk)
                created_markers += 1

        if not dry_run:
            db.commit()
        else:
            db.rollback()

    return created_scd, created_clusters, created_markers


def main() -> None:
    args = _parse_args()
    if args.reset:
        ds_deleted, sc_deleted = _reset_demo()
        print(f"Reset: deleted demo Datasets = {ds_deleted}")
        print(f"Reset: deleted SingleCellDataset rows = {sc_deleted}")

    created_scd, created_clusters, created_markers = _seed(args.size, args.seed, args.dry_run)
    print(f"Seed complete (size={args.size}, dry_run={args.dry_run})")
    print(f"SingleCellDataset created: {created_scd}")
    print(f"CellCluster created: {created_clusters}")
    print(f"CellTypeMarker created: {created_markers}")


if __name__ == "__main__":
    main()


