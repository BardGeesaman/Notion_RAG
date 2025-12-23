"""Multi-omics latent factor analysis service (MOFA orchestration)."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple
from uuid import UUID

import pandas as pd

from amprenta_rag.database.models import (
    FactorLoading,
    FactorScore,
    Feature,
    LatentFactor,
    MultiOmicsExperiment,
)
from amprenta_rag.multi_omics.data_extractor import align_samples, extract_feature_matrix
from amprenta_rag.multi_omics.mofa_runner import run_mofa
from amprenta_rag.multi_omics.result_parser import parse_mofa_output


def _dataset_entries(exp: MultiOmicsExperiment, db) -> List[Tuple[str, UUID]]:
    """Resolve experiment.dataset_ids JSON into (omics_type, dataset_id) tuples."""
    ds_ids = exp.dataset_ids or []
    out: List[Tuple[str, UUID]] = []

    def _as_uuid(x: Any) -> UUID:
        if isinstance(x, UUID):
            return x
        return UUID(str(x))

    if isinstance(ds_ids, dict):
        for k, v in ds_ids.items():
            out.append((str(k), _as_uuid(v)))
        return out

    if isinstance(ds_ids, list):
        for item in ds_ids:
            if isinstance(item, dict):
                did = item.get("dataset_id") or item.get("id")
                if not did:
                    continue
                ot = item.get("omics_type") or item.get("view") or item.get("omics")
                did_u = _as_uuid(did)
                if not ot:
                    ot = _infer_omics_type(did_u, db)
                out.append((str(ot), did_u))
            else:
                did_u = _as_uuid(item)
                out.append((_infer_omics_type(did_u, db), did_u))
        return out

    raise ValueError("MultiOmicsExperiment.dataset_ids must be list or dict")


def _infer_omics_type(dataset_id: UUID, db) -> str:
    from amprenta_rag.database.models import Dataset as DatasetModel

    ds = db.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    if not ds:
        return "omics"
    ot = getattr(ds, "omics_type", None) or "omics"
    return str(ot)


def _map_features_to_ids(feature_names: List[str], db) -> Dict[str, UUID]:
    names = [n for n in feature_names if n]
    if not names:
        return {}
    rows = db.query(Feature).filter(Feature.name.in_(names)).all()
    out: Dict[str, UUID] = {r.name: r.id for r in rows if r and r.name}
    return out


def _variance_per_factor(
    variance: Dict[str, Any], views: List[str], n_factors: int
) -> Dict[int, Dict[str, Any]]:
    """Best-effort parse of variance explained into per-factor dicts {factor_i: {view: value}}."""
    out: Dict[int, Dict[str, Any]] = {i: {} for i in range(n_factors)}
    for v in variance.values():
        if not isinstance(v, list) or not v:
            continue
        if not all(isinstance(row, list) for row in v):
            continue
        if len(v) != len(views):
            continue
        if not all(len(row) == n_factors for row in v):
            continue
        for fi in range(n_factors):
            out[fi] = {views[vi]: v[vi][fi] for vi in range(len(views))}
        return out
    return out


def run_experiment(experiment_id: UUID, db) -> List[LatentFactor]:
    """Run MOFA for a MultiOmicsExperiment and persist factors/loadings/scores."""
    exp = db.query(MultiOmicsExperiment).filter(MultiOmicsExperiment.id == experiment_id).first()
    if not exp:
        raise ValueError("MultiOmicsExperiment not found")

    log_lines: List[str] = []
    try:
        exp.status = "running"
        exp.processing_log = None
        db.add(exp)
        db.commit()

        # Resolve datasets
        entries = _dataset_entries(exp, db)
        matrices: Dict[str, pd.DataFrame] = {}
        for ot, did in entries:
            omics_type = ot if ot and ot != "omics" else _infer_omics_type(did, db)
            key = str(omics_type)
            if key in matrices:
                key = f"{key}:{str(did)[:8]}"
            mat = extract_feature_matrix(did, db)
            matrices[key] = mat
            log_lines.append(f"Loaded matrix {key} shape={mat.shape}")

        aligned = align_samples(matrices, exp.sample_mapping or [])
        for k, m in aligned.items():
            log_lines.append(f"Aligned {k} shape={m.shape}")

        out_dir = Path("data") / "multi_omics" / str(experiment_id)
        out_dir.mkdir(parents=True, exist_ok=True)
        h5_path = run_mofa(
            aligned,
            n_factors=int(exp.n_factors or 10),
            convergence_mode=str(exp.convergence_mode or "fast"),
            output_dir=str(out_dir),
        )
        log_lines.append(f"MOFA output: {h5_path}")

        parsed = parse_mofa_output(h5_path)
        factors = parsed.get("factors") or []
        scores_df: pd.DataFrame = parsed.get("scores")
        loadings: Dict[str, pd.DataFrame] = parsed.get("loadings") or {}
        variance: Dict[str, Any] = parsed.get("variance") or {}

        # Clear existing results (manual deletes to be safe with bulk operations)
        factor_ids = [r[0] for r in db.query(LatentFactor.id).filter(LatentFactor.experiment_id == experiment_id).all()]
        if factor_ids:
            db.query(FactorLoading).filter(FactorLoading.factor_id.in_(factor_ids)).delete(synchronize_session=False)
            db.query(FactorScore).filter(FactorScore.factor_id.in_(factor_ids)).delete(synchronize_session=False)
            db.query(LatentFactor).filter(LatentFactor.experiment_id == experiment_id).delete(synchronize_session=False)
            db.commit()

        created: List[LatentFactor] = []

        # Preload feature mappings per view (by Feature.name)
        feature_maps: Dict[str, Dict[str, UUID]] = {}
        for view, wdf in loadings.items():
            feature_maps[view] = _map_features_to_ids(list(map(str, wdf.index.tolist())), db)
            log_lines.append(f"Mapped features for {view}: {len(feature_maps[view])}/{wdf.shape[0]}")

        # Build per-factor variance dict if variance contains view->list
        view_names = list(loadings.keys())
        variance_by_factor = _variance_per_factor(variance, views=view_names, n_factors=len(factors))

        # Create factors + persist loadings and scores
        for i in factors:
            lf = LatentFactor(
                experiment_id=experiment_id,
                factor_index=int(i),
                variance_explained=variance_by_factor.get(int(i)) or None,
                description=None,
            )
            db.add(lf)
            db.flush()  # get lf.id

            # Scores: scores_df index = samples, columns = factors
            if scores_df is not None and int(i) in list(scores_df.columns):
                for sid, val in scores_df[int(i)].items():
                    try:
                        sc = float(val)
                    except Exception:
                        continue
                    db.add(FactorScore(factor_id=lf.id, sample_id=str(sid), score=sc))

            # Loadings per view: wdf index=features, columns=factors
            for view, wdf in loadings.items():
                fmap = feature_maps.get(view) or {}
                if int(i) not in list(wdf.columns):
                    continue
                col = wdf[int(i)]
                for feat_name, val in col.items():
                    fid = fmap.get(str(feat_name))
                    if not fid:
                        continue
                    try:
                        lv = float(val)
                    except Exception:
                        continue
                    db.add(
                        FactorLoading(
                            factor_id=lf.id,
                            feature_id=fid,
                            loading=lv,
                            abs_loading=abs(lv),
                            omics_type=str(view),
                        )
                    )

            created.append(lf)

        exp.status = "completed"
        exp.processed_at = datetime.utcnow()
        exp.processing_log = "\n".join(log_lines)[-20000:] if log_lines else None
        db.add(exp)
        db.commit()

        return created
    except Exception as e:  # noqa: BLE001
        exp.status = "failed"
        exp.processed_at = datetime.utcnow()
        msg = "\n".join(log_lines + [f"ERROR: {e}"])[-20000:]
        exp.processing_log = msg
        db.add(exp)
        db.commit()
        raise


__all__ = ["run_experiment"]


