"""LINCS ingestion orchestration (download -> parse -> DB insert)."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterator

from amprenta_rag.connectivity.lincs_fetcher import download_lincs_level5
from amprenta_rag.connectivity.lincs_parser import LINCSSignatureData, parse_gct, parse_gctx
from amprenta_rag.database.models import LINCSGene, LINCSSignature
from amprenta_rag.database.session import db_session


def _iter_signatures(path: str) -> Iterator[LINCSSignatureData]:
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".gct":
        yield from parse_gct(path)
    elif suf == ".gctx":
        yield from parse_gctx(path)
    else:
        raise ValueError(f"Unsupported LINCS file type: {p.suffix}")


def ingest_lincs_data(max_compounds: int = 1000) -> int:
    """Download and ingest LINCS signatures into DB.

    For MVP, ingests up to `max_compounds` signatures (columns) from the parsed file.
    """
    base = Path(os.environ.get("LINCS_STORE_BASE", "data/lincs"))
    base.mkdir(parents=True, exist_ok=True)

    path = download_lincs_level5(str(base), max_compounds=max_compounds)
    ingested = 0
    seen_entrez: set[int] = set()

    with db_session() as db:
        for sig in _iter_signatures(path):
            # Limit
            if ingested >= int(max_compounds):
                break

            # Upsert-ish for signature by sig_id
            existing = db.query(LINCSSignature).filter_by(sig_id=sig.sig_id).first()
            if existing is None:
                srow = LINCSSignature(
                    sig_id=sig.sig_id,
                    pert_iname=sig.pert_iname,
                    pert_id=sig.pert_id,
                    pert_type=None,
                    cell_id=sig.cell_id,
                    pert_time=None,
                    pert_dose=None,
                    gene_expression={str(k): float(v) for k, v in (sig.gene_expression or {}).items()},
                    tas=None,
                    distil_id=None,
                )
                db.add(srow)
            else:
                existing.gene_expression = {str(k): float(v) for k, v in (sig.gene_expression or {}).items()}
                existing.pert_iname = sig.pert_iname
                existing.pert_id = sig.pert_id
                existing.cell_id = sig.cell_id

            # Create missing genes as placeholders
            for entrez in (sig.gene_expression or {}).keys():
                if entrez in seen_entrez:
                    continue
                seen_entrez.add(entrez)
                if db.query(LINCSGene).filter_by(entrez_id=int(entrez)).first() is None:
                    db.add(LINCSGene(entrez_id=int(entrez), gene_symbol=None, gene_title=None, is_landmark=False))

            ingested += 1

    return ingested


__all__ = ["ingest_lincs_data"]


