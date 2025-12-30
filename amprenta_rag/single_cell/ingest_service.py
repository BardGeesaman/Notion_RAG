"""Single-cell ingestion service (h5ad -> preprocess -> markers -> DB tables)."""

from __future__ import annotations

import threading
from datetime import datetime, timezone
from pathlib import Path
from uuid import UUID

from amprenta_rag.database.models import (
    CellAnnotation,
    CellCluster,
    CellTypeMarker,
    Dataset,
    SingleCellDataset,
)
from amprenta_rag.database.session import db_session
from amprenta_rag.single_cell.gene_mapper import map_genes_to_features
from amprenta_rag.single_cell.h5ad_parser import extract_metadata, load_h5ad, validate_h5ad
from amprenta_rag.single_cell.marker_discovery import find_markers
from amprenta_rag.single_cell.scanpy_pipeline import run_preprocessing


def ingest_h5ad(h5ad_path: str, dataset_id: UUID | None = None) -> SingleCellDataset:
    """Create SingleCellDataset row and start background processing."""
    p = Path(h5ad_path)
    if not p.exists():
        raise FileNotFoundError(str(p))

    with db_session() as db:
        ds_id = dataset_id
        if ds_id is None:
            # Minimal Dataset row so we can FK; omics_type kept generic.
            ds = Dataset(name=p.stem, omics_type="single_cell", description="Imported from h5ad")
            db.add(ds)
            db.commit()
            db.refresh(ds)
            ds_id = ds.id

        scd = SingleCellDataset(
            dataset_id=ds_id,
            h5ad_path=str(p),
            file_size_bytes=p.stat().st_size,
            processing_status="pending",
            processing_log=None,
            ingested_at=datetime.now(timezone.utc),
        )
        db.add(scd)
        db.commit()
        db.refresh(scd)

        scd_id = scd.id

    thread = threading.Thread(target=_process_single_cell_async, args=(scd_id,), daemon=True)
    thread.start()

    with db_session() as db:
        out = db.query(SingleCellDataset).filter_by(id=scd_id).first()
        if not out:
            raise RuntimeError("Failed to load created SingleCellDataset")
        return out


def recluster_dataset(single_cell_dataset_id: UUID, *, resolution: float = 0.5) -> None:
    """Trigger a re-clustering run in the background (MVP: full reprocess)."""
    with db_session() as db:
        scd = db.query(SingleCellDataset).filter_by(id=single_cell_dataset_id).first()
        if not scd:
            raise ValueError("SingleCellDataset not found")
        scd.clustering_resolution = float(resolution)
        scd.processing_status = "pending"
        scd.processing_log = None
        db.add(scd)

    # Submit job via Celery or fallback to threading
    import os
    if os.environ.get("USE_CELERY", "true").lower() == "true":
        from amprenta_rag.jobs.tasks.single_cell import process_single_cell
        process_single_cell.delay(str(single_cell_dataset_id))
    else:
        # Fallback to threading for gradual rollout
        thread = threading.Thread(target=_process_single_cell_async, args=(single_cell_dataset_id,), daemon=True)
        thread.start()


def _process_single_cell_async(single_cell_dataset_id: UUID) -> None:
    """Background processing: preprocess + markers + populate tables."""
    try:
        with db_session() as db:
            scd = db.query(SingleCellDataset).filter_by(id=single_cell_dataset_id).first()
            if not scd:
                return
            scd.processing_status = "running"
            scd.processing_log = None
            db.add(scd)

        adata = load_h5ad(scd.h5ad_path)
        validate_h5ad(adata)
        meta = extract_metadata(adata)

        resolution = float(getattr(scd, "clustering_resolution", None) or 0.5)
        processed = run_preprocessing(adata, resolution=resolution)
        markers = find_markers(processed, groupby="leiden", top_n=50)

        # Build per-cell annotations
        obs = processed.obs
        umap = processed.obsm.get("X_umap")

        with db_session() as db:
            scd = db.query(SingleCellDataset).filter_by(id=single_cell_dataset_id).first()
            if not scd:
                return

            scd.n_cells = int(meta.get("n_cells") or processed.n_obs)
            scd.n_genes = int(meta.get("n_genes") or processed.n_vars)
            scd.normalization_method = "normalize_total+log1p"
            scd.hvg_count = int(processed.n_vars)
            scd.clustering_method = "leiden"
            scd.clustering_resolution = float(scd.clustering_resolution or resolution)
            scd.has_pca = True
            scd.has_umap = umap is not None

            # Clear existing derived tables (idempotent re-run)
            db.query(CellAnnotation).filter_by(single_cell_dataset_id=single_cell_dataset_id).delete()
            db.query(CellCluster).filter_by(single_cell_dataset_id=single_cell_dataset_id).delete()
            db.query(CellTypeMarker).filter_by(single_cell_dataset_id=single_cell_dataset_id).delete()

            # Insert annotations
            for i, barcode in enumerate(list(processed.obs_names)):
                cluster = None
                try:
                    cluster = int(obs.loc[barcode, "leiden"])
                except Exception:
                    cluster = None
                cell_type = obs.loc[barcode].get("cell_type") if hasattr(obs.loc[barcode], "get") else None
                u1 = float(umap[i, 0]) if umap is not None else None
                u2 = float(umap[i, 1]) if umap is not None else None

                db.add(
                    CellAnnotation(
                        single_cell_dataset_id=single_cell_dataset_id,
                        barcode=str(barcode),
                        cluster_id=cluster,
                        cell_type=str(cell_type) if cell_type else None,
                        umap_1=u1,
                        umap_2=u2,
                        n_genes_detected=None,
                        total_counts=None,
                        pct_mitochondrial=None,
                    )
                )

            # Cluster summaries
            if "leiden" in obs:
                counts = obs["leiden"].value_counts()
                for cid, n_cells in counts.items():
                    try:
                        cluster_id = int(cid)
                    except Exception:
                        continue
                    db.add(
                        CellCluster(
                            single_cell_dataset_id=single_cell_dataset_id,
                            cluster_id=cluster_id,
                            cell_type=None,
                            n_cells=int(n_cells),
                            marker_feature_ids=None,
                            description=None,
                        )
                    )

            # Marker ingestion
            gene_symbols = [str(x) for x in markers.get("gene_symbol", []).tolist()] if not markers.empty else []
            gene_to_feature = map_genes_to_features(gene_symbols, db)

            for _, row in markers.iterrows():
                gene = str(row.get("gene_symbol"))
                cid = row.get("cluster_id")
                try:
                    cluster_id = int(cid)
                except Exception:
                    continue
                fid = gene_to_feature.get(gene)
                db.add(
                    CellTypeMarker(
                        single_cell_dataset_id=single_cell_dataset_id,
                        cluster_id=cluster_id,
                        feature_id=fid,
                        gene_symbol=gene,
                        log2_fold_change=float(row.get("log2_fold_change")) if row.get("log2_fold_change") is not None else None,
                        pval=float(row.get("pval")) if row.get("pval") is not None else None,
                        pval_adj=float(row.get("pval_adj")) if row.get("pval_adj") is not None else None,
                        pct_in_cluster=float(row.get("pct_in_cluster")) if row.get("pct_in_cluster") is not None else None,
                        pct_out_cluster=float(row.get("pct_out_cluster")) if row.get("pct_out_cluster") is not None else None,
                    )
                )

            scd.processing_status = "completed"
            scd.processed_at = datetime.now(timezone.utc)
            db.add(scd)

    except Exception as e:  # noqa: BLE001
        with db_session() as db:
            scd = db.query(SingleCellDataset).filter_by(id=single_cell_dataset_id).first()
            if scd:
                scd.processing_status = "failed"
                scd.processing_log = str(e)[:4000]
                scd.processed_at = datetime.now(timezone.utc)
                db.add(scd)


__all__ = ["ingest_h5ad"]


