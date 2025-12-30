"""Celery tasks for single-cell processing."""

from datetime import datetime, timezone
from uuid import UUID

from amprenta_rag.jobs.celery_app import celery_app


@celery_app.task(bind=True, max_retries=2, default_retry_delay=60, queue='default')
def process_single_cell(self, dataset_id: str) -> dict:
    """Process single-cell dataset (QC, normalize, cluster)."""
    from amprenta_rag.database.session import db_session
    from amprenta_rag.database.models import (
        CellAnnotation,
        CellCluster,
        CellTypeMarker,
        SingleCellDataset,
    )
    from amprenta_rag.single_cell.gene_mapper import map_genes_to_features
    from amprenta_rag.single_cell.h5ad_parser import extract_metadata, load_h5ad, validate_h5ad
    from amprenta_rag.single_cell.marker_discovery import find_markers
    from amprenta_rag.single_cell.scanpy_pipeline import run_preprocessing
    
    dataset_uuid = UUID(dataset_id)
    
    try:
        # Mark running
        with db_session() as db:
            scd = db.query(SingleCellDataset).filter_by(id=dataset_uuid).first()
            if not scd:
                return {"status": "failed", "error": "Dataset not found", "dataset_id": dataset_id}
            scd.processing_status = "running"
            scd.processing_log = None
            db.add(scd)
            db.commit()

        # Load and validate h5ad
        adata = load_h5ad(scd.h5ad_path)
        validate_h5ad(adata)
        meta = extract_metadata(adata)

        # Run preprocessing pipeline
        resolution = float(getattr(scd, "clustering_resolution", None) or 0.5)
        processed = run_preprocessing(adata, resolution=resolution)
        markers = find_markers(processed, groupby="leiden", top_n=50)

        # Build per-cell annotations
        obs = processed.obs
        umap = processed.obsm.get("X_umap")

        with db_session() as db:
            scd = db.query(SingleCellDataset).filter_by(id=dataset_uuid).first()
            if not scd:
                return {"status": "failed", "error": "Dataset not found", "dataset_id": dataset_id}

            # Update dataset metadata
            scd.n_cells = int(meta.get("n_cells") or processed.n_obs)
            scd.n_genes = int(meta.get("n_genes") or processed.n_vars)
            scd.normalization_method = "normalize_total+log1p"
            scd.hvg_count = int(processed.n_vars)
            scd.clustering_method = "leiden"
            scd.clustering_resolution = float(scd.clustering_resolution or resolution)
            scd.has_pca = True
            scd.has_umap = umap is not None

            # Clear existing derived tables (idempotent re-run)
            db.query(CellAnnotation).filter_by(single_cell_dataset_id=dataset_uuid).delete()
            db.query(CellCluster).filter_by(single_cell_dataset_id=dataset_uuid).delete()
            db.query(CellTypeMarker).filter_by(single_cell_dataset_id=dataset_uuid).delete()

            # Insert cell annotations
            for i, barcode in enumerate(obs.index):
                leiden_cluster = obs.iloc[i].get("leiden")
                try:
                    cluster_id = int(leiden_cluster)
                except Exception:
                    cluster_id = 0

                umap_x = umap[i, 0] if umap is not None else None
                umap_y = umap[i, 1] if umap is not None else None

                db.add(
                    CellAnnotation(
                        single_cell_dataset_id=dataset_uuid,
                        cell_barcode=str(barcode),
                        cluster_id=cluster_id,
                        umap_x=float(umap_x) if umap_x is not None else None,
                        umap_y=float(umap_y) if umap_y is not None else None,
                    )
                )

            # Insert cluster summaries
            for cluster_id in obs["leiden"].unique():
                try:
                    cid = int(cluster_id)
                except Exception:
                    continue
                
                cluster_cells = (obs["leiden"] == cluster_id).sum()
                db.add(
                    CellCluster(
                        single_cell_dataset_id=dataset_uuid,
                        cluster_id=cid,
                        cluster_label=f"Cluster {cid}",
                        cell_count=int(cluster_cells),
                        cluster_color=None,
                    )
                )

            # Insert markers
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
                        single_cell_dataset_id=dataset_uuid,
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

            # Mark complete
            scd.processing_status = "completed"
            scd.processed_at = datetime.now(timezone.utc)
            db.add(scd)
            db.commit()

        return {
            "status": "completed",
            "dataset_id": dataset_id,
            "n_cells": scd.n_cells,
            "n_genes": scd.n_genes,
            "n_clusters": len(obs["leiden"].unique())
        }

    except Exception as exc:
        # Update dataset status to failed
        try:
            with db_session() as db:
                scd = db.query(SingleCellDataset).filter_by(id=dataset_uuid).first()
                if scd is not None:
                    scd.processing_status = "failed"
                    scd.processing_log = str(exc)[:4000]
                    scd.processed_at = datetime.now(timezone.utc)
                    db.add(scd)
                    db.commit()
        except Exception:
            pass  # Don't let DB errors prevent retry logic
        
        # Retry with exponential backoff
        if self.request.retries >= self.max_retries:
            return {"status": "failed", "error": str(exc), "dataset_id": dataset_id}
        
        self.retry(exc=exc, countdown=60 * (2 ** self.request.retries))
