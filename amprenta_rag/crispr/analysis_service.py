"""CRISPR screen analysis service (MAGeCK)."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional
from uuid import UUID

from amprenta_rag.crispr.count_parser import parse_count_matrix
from amprenta_rag.crispr.gene_mapper import map_genes_to_features
from amprenta_rag.crispr.mageck_runner import run_mageck_test
from amprenta_rag.crispr.result_parser import parse_gene_summary
from amprenta_rag.database.models import CRISPRResult, CRISPRScreen, Dataset


def _count_file_for_screen(screen: CRISPRScreen) -> str:
    ds = getattr(screen, "dataset", None)
    if ds is None:
        raise ValueError("CRISPRScreen.dataset not loaded")
    fp = (ds.file_paths or [None])[0]
    if not fp:
        raise ValueError(f"Dataset {ds.id} has no file_paths")
    return str(fp)


def _ensure_samples_present(count_path: str, control: str, treatment: str) -> None:
    df = parse_count_matrix(count_path)
    cols = set(df.columns)
    if control not in cols:
        raise ValueError(f"Control sample '{control}' not found in count matrix columns")
    if treatment not in cols:
        raise ValueError(f"Treatment sample '{treatment}' not found in count matrix columns")


def run_screen_analysis(
    screen_id: UUID,
    db,
    *,
    method: str = "test",
) -> List[CRISPRResult]:
    """Run CRISPR screen analysis and persist CRISPRResult rows.

    - Loads CRISPRScreen + linked Dataset.
    - Uses Dataset.file_paths[0] as the MAGeCK count matrix file.
    - Runs MAGeCK, parses gene_summary.
    - Maps gene_symbol -> Feature.id (feature_type='gene').
    - Clears existing CRISPRResult rows for the screen and inserts new results.
    - Sets is_hit if fdr < 0.05.
    - Updates screen.status to completed or failed.
    """
    if method != "test":
        raise ValueError("Only method='test' is supported for MAGeCK runner MVP")

    screen = db.query(CRISPRScreen).filter(CRISPRScreen.id == screen_id).first()
    if not screen:
        raise ValueError("CRISPRScreen not found")

    # Ensure dataset loaded
    screen.dataset = db.query(Dataset).filter(Dataset.id == screen.dataset_id).first()
    if not screen.dataset:
        raise ValueError("Linked Dataset not found")

    control = (screen.control_label or "").strip()
    treatment = (screen.treatment_label or "").strip()
    if not control or not treatment:
        raise ValueError("CRISPRScreen.control_label and treatment_label are required")

    try:
        screen.status = "running"
        db.add(screen)
        db.commit()

        count_file = _count_file_for_screen(screen)
        _ensure_samples_present(count_file, control=control, treatment=treatment)

        out_dir = Path("data") / "crispr" / str(screen_id)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = out_dir / "mageck"

        gene_summary_path = run_mageck_test(
            count_file=count_file,
            control=control,
            treatment=treatment,
            output_prefix=str(out_prefix),
        )

        parsed = parse_gene_summary(gene_summary_path)
        genes = [str(r.get("gene")) for r in parsed if r.get("gene")]
        gene_to_feature = map_genes_to_features(genes, db)

        # Clear existing results and re-insert
        db.query(CRISPRResult).filter(CRISPRResult.screen_id == screen_id).delete()

        # Rank by FDR ascending, fall back to large number
        def _fdr_key(r) -> float:
            v = r.get("fdr")
            try:
                return float(v) if v is not None else 1e9
            except Exception:
                return 1e9

        parsed_sorted = sorted(parsed, key=_fdr_key)

        out: List[CRISPRResult] = []
        for i, r in enumerate(parsed_sorted, start=1):
            gene = r.get("gene")
            if gene is None:
                continue

            neg_p = r.get("neg_p")
            pos_p = r.get("pos_p")
            p_value: Optional[float]
            try:
                pv = [float(x) for x in (neg_p, pos_p) if x is not None]
                p_value = min(pv) if pv else None
            except Exception:
                p_value = None

            fdr = r.get("fdr")
            try:
                fdr_f = float(fdr) if fdr is not None else None
            except Exception:
                fdr_f = None

            feat_id = gene_to_feature.get(str(gene))

            row = CRISPRResult(
                screen_id=screen_id,
                gene_symbol=str(gene),
                feature_id=feat_id,
                beta_score=None,
                p_value=p_value,
                fdr=fdr_f,
                neg_lfc=r.get("neg_lfc"),
                pos_lfc=r.get("pos_lfc"),
                rank=i,
                is_hit=bool(fdr_f is not None and fdr_f < 0.05),
            )
            db.add(row)
            out.append(row)

        screen.status = "completed"
        db.add(screen)
        db.commit()

        return out
    except Exception:
        screen.status = "failed"
        db.add(screen)
        db.commit()
        raise


__all__ = ["run_screen_analysis"]


