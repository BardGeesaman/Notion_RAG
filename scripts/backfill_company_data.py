from __future__ import annotations

import argparse
from typing import Iterable, List, Optional

from sqlalchemy import text

from amprenta_rag.database.models import Company
from amprenta_rag.database.session import db_session


TENANT_TABLES = [
    "datasets",
    "signatures",
    "experiments",
    "programs",
    "compounds",
    "hts_campaigns",
    "hts_plates",
    "variant_sets",
    "single_cell_datasets",
    "crispr_screens",
    "multi_omics_experiments",
    "ml_models",
    "docking_runs",
    "batch_corrections",
]


def _table_columns(db, table: str) -> List[str]:
    rows = db.execute(
        text(
            """
SELECT column_name
FROM information_schema.columns
WHERE table_schema = 'public' AND table_name = :t
"""
        ),
        {"t": table},
    ).fetchall()
    return [r[0] for r in rows]


def _get_default_company_id(db, subdomain: str) -> str:
    comp = db.query(Company).filter(Company.subdomain == subdomain).first()
    if not comp:
        raise RuntimeError(f"Default company with subdomain='{subdomain}' not found. Run scripts/seed_default_company.py first.")
    return str(comp.id)


def _backfill_table(db, table: str, default_company_id: str, *, dry_run: bool) -> None:
    cols = set(_table_columns(db, table))
    if "company_id" not in cols:
        return

    # Strategy priority:
    # 1) created_by_id / created_by / uploaded_by -> users.company_id
    # 2) dataset_id -> datasets.company_id
    # 3) fallback -> default company
    creator_col = None
    for cand in ("created_by_id", "created_by", "uploaded_by"):
        if cand in cols:
            creator_col = cand
            break

    if creator_col:
        stmt = text(
            f"""
UPDATE "{table}" t
SET company_id = COALESCE(u.company_id, :default_company_id::uuid)
FROM users u
WHERE t.company_id IS NULL
  AND t.{creator_col} IS NOT NULL
  AND u.id = t.{creator_col};
"""
        )
        if not dry_run:
            db.execute(stmt, {"default_company_id": default_company_id})

    elif "dataset_id" in cols:
        stmt = text(
            f"""
UPDATE "{table}" t
SET company_id = d.company_id
FROM datasets d
WHERE t.company_id IS NULL
  AND t.dataset_id IS NOT NULL
  AND d.id = t.dataset_id
  AND d.company_id IS NOT NULL;
"""
        )
        if not dry_run:
            db.execute(stmt)

    # Final fallback: assign default company for any remaining NULL company_id
    stmt2 = text(
        f"""
UPDATE "{table}"
SET company_id = :default_company_id::uuid
WHERE company_id IS NULL;
"""
    )
    if not dry_run:
        db.execute(stmt2, {"default_company_id": default_company_id})


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Backfill company_id for tenant tables.")
    ap.add_argument("--default-subdomain", type=str, default="default")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(list(argv) if argv is not None else None)

    with db_session() as db:
        default_company_id = _get_default_company_id(db, args.default_subdomain)
        for t in TENANT_TABLES:
            try:
                _backfill_table(db, t, default_company_id, dry_run=bool(args.dry_run))
            except Exception as e:  # noqa: BLE001
                print(f"[WARN] backfill failed for table={t}: {e}")
        if not args.dry_run:
            db.commit()

    print(f"Backfill complete. default_company_id={default_company_id} dry_run={args.dry_run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


