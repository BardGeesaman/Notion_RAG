"""Load HPO phenotype->gene associations into the database.

Source:
  http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt

This script downloads the file (if needed), parses it, and upserts into Postgres.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import requests
from sqlalchemy.dialects.postgresql import insert

from amprenta_rag.database.models import PhenotypeGeneAssociation
from amprenta_rag.database.session import db_session
from amprenta_rag.phenotype.hpo_parser import parse_genes_to_phenotype


DEFAULT_URL = "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"


def download(url: str, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    out_path.write_bytes(r.content)
    return out_path


def load_records(path: Path, limit: int | None = None) -> int:
    records = parse_genes_to_phenotype(path)
    if limit is not None:
        records = records[:limit]

    rows: List[Dict[str, object]] = []
    for r in records:
        rows.append(
            {
                "id": None,  # let ORM default if used; for SQL insert we omit it
                "hpo_id": r.hpo_id,
                "hpo_name": r.hpo_name,
                "gene_symbol": r.gene_symbol,
                "ncbi_gene_id": r.ncbi_gene_id,
                "disease_id": r.disease_id,
                "frequency": r.frequency,
                "source": r.source,
            }
        )

    with db_session() as db:
        created = 0
        # Chunked bulk insert with ON CONFLICT DO NOTHING (Postgres).
        try:
            tbl = PhenotypeGeneAssociation.__table__
            chunk_size = 5000
            for i in range(0, len(rows), chunk_size):
                chunk = rows[i : i + chunk_size]
                # Remove "id" to let default UUID generator handle it (db-side / ORM).
                chunk2 = [{k: v for k, v in r.items() if k != "id"} for r in chunk]
                stmt = insert(tbl).values(chunk2).on_conflict_do_nothing(constraint="uq_hpo_gene_disease")
                res = db.execute(stmt)
                # rowcount is best-effort with ON CONFLICT DO NOTHING; can be -1 depending on driver.
                if res.rowcount and res.rowcount > 0:
                    created += int(res.rowcount)
            db.commit()
            return created
        except Exception:
            # Fallback: ORM inserts with per-row merge-like behavior (slower but portable).
            for r in records:
                obj = PhenotypeGeneAssociation(
                    hpo_id=r.hpo_id,
                    hpo_name=r.hpo_name,
                    gene_symbol=r.gene_symbol,
                    ncbi_gene_id=r.ncbi_gene_id,
                    disease_id=r.disease_id,
                    frequency=r.frequency,
                    source=r.source,
                )
                try:
                    db.add(obj)
                    db.flush()
                    created += 1
                except Exception:
                    db.rollback()
            db.commit()
            return created


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--out", default="data/hpo/genes_to_phenotype.txt")
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    out_path = Path(args.out)
    if not args.no_download and not out_path.exists():
        print(f"Downloading {args.url} -> {out_path}")
        download(args.url, out_path)

    print(f"Loading {out_path} ...")
    created = load_records(out_path, limit=args.limit)
    print(f"Inserted {created} phenotype-gene associations.")


if __name__ == "__main__":
    main()


